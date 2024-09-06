use std::fmt::Debug;
use std::marker::PhantomData;

use ark_ff::{BigInteger, PrimeField};
use p3_challenger::{DuplexChallenger, HashChallenger, SerializingChallenger32};
use p3_circle::CirclePcs;
use p3_commit::ExtensionMmcs;
use p3_dft::Radix2DitParallel;
use p3_field::extension::BinomialExtensionField;
use p3_field::{AbstractField, Field};
use p3_fri::{FriConfig, TwoAdicFriPcs};
use p3_goldilocks::{DiffusionMatrixGoldilocks, Goldilocks};
use p3_keccak::Keccak256Hash;
use p3_merkle_tree::FieldMerkleTreeMmcs;
use p3_poseidon2::{Poseidon2, Poseidon2ExternalMatrixGeneral};
use p3_poseidon2_air::{generate_trace, FieldType, Poseidon2Air};
use p3_symmetric::{
    CompressionFunctionFromHasher, PaddingFreeSponge, SerializingHasher32, TruncatedPermutation,
};
use p3_uni_stark::{prove, verify, StarkConfig};
use rand::{random, thread_rng};
use tracing_forest::util::LevelFilter;
use tracing_forest::ForestLayer;
use tracing_subscriber::layer::SubscriberExt;
use tracing_subscriber::util::SubscriberInitExt;
use tracing_subscriber::{EnvFilter, Registry};
use zkhash::fields::goldilocks::FpGoldiLocks;
use zkhash::poseidon::poseidon_instance_goldilocks::RC16;

const WIDTH: usize = 16;
const SBOX_DEGREE: u64 = 3;
const SBOX_REGISTERS: usize = 1;
const HALF_FULL_ROUNDS: usize = 4;
const PARTIAL_ROUNDS: usize = 20;

const NUM_HASHES: usize = 1 << 4;

fn goldilocks_from_ark_ff(input: FpGoldiLocks) -> Goldilocks {
    let as_bigint = input.into_bigint();
    let mut as_bytes = as_bigint.to_bytes_le();
    as_bytes.resize(8, 0);
    let as_u64 = u64::from_le_bytes(as_bytes[0..8].try_into().unwrap());
    Goldilocks::from_wrapped_u64(as_u64)
}

fn main() {
    let env_filter = EnvFilter::builder()
        .with_default_directive(LevelFilter::INFO.into())
        .from_env_lossy();

    Registry::default()
        .with(env_filter)
        .with(ForestLayer::default())
        .init();

    type Val = Goldilocks;
    type Challenge = BinomialExtensionField<Val, 2>;
    type Perm = Poseidon2<
        Goldilocks,
        Poseidon2ExternalMatrixGeneral,
        DiffusionMatrixGoldilocks,
        WIDTH,
        SBOX_DEGREE,
    >;

    let round_constants: Vec<[Val; WIDTH]> = RC16
        .iter()
        .map(|vec| {
            vec.iter()
                .cloned()
                .map(goldilocks_from_ark_ff)
                .collect::<Vec<_>>()
                .try_into()
                .unwrap()
        })
        .collect();

    let external_linear_layer = Poseidon2ExternalMatrixGeneral::default();
    let internal_linear_layer = DiffusionMatrixGoldilocks::default();

    let mut external_constants = round_constants.clone()[0..4].to_vec();
    external_constants.extend(round_constants.clone()[26..30].to_vec());

    let mut internal_constants = round_constants.clone()[4..26]
        .iter()
        .map(|x| x[0])
        .collect();

    let perm = Perm::new(
        8,
        external_constants,
        external_linear_layer,
        22,
        internal_constants,
        internal_linear_layer,
    );

    type MyHash = PaddingFreeSponge<Perm, WIDTH, 4, 4>;
    let hash = MyHash::new(perm.clone());

    type MyCompress = TruncatedPermutation<Perm, 2, 4, WIDTH>;
    let compress = MyCompress::new(perm.clone());

    type ValMmcs = FieldMerkleTreeMmcs<
        <Val as Field>::Packing,
        <Val as Field>::Packing,
        MyHash,
        MyCompress,
        4,
    >;
    let val_mmcs = ValMmcs::new(hash, compress);

    type ChallengeMmcs = ExtensionMmcs<Val, Challenge, ValMmcs>;
    let challenge_mmcs = ChallengeMmcs::new(val_mmcs.clone());

    type Dft = Radix2DitParallel;
    let dft = Dft {};

    type Challenger = DuplexChallenger<Val, Perm, WIDTH, 4>;

    //     let external_linear_layer = Poseidon2ExternalMatrixGeneral::default();
    //     let internal_linear_layer = DiffusionMatrixGoldilocks::default();
    let air: Poseidon2Air<Val, WIDTH> = Poseidon2Air::new();
    let mut input = core::array::from_fn(|j| Goldilocks::from_canonical_u64(0));
    let trace = generate_trace::<Val, WIDTH>(&mut input, round_constants.clone(), FieldType::GL64);

    let fri_config = FriConfig {
        log_blowup: 1,
        num_queries: 100,
        proof_of_work_bits: 16,
        mmcs: challenge_mmcs,
    };
    type Pcs = TwoAdicFriPcs<Val, Dft, ValMmcs, ChallengeMmcs>;

    // dbg!(log2_ceil_usize(trace.height()));
    // dbg!(get_log_quotient_degree::<Val, FibonacciAir>(
    //     &FibonacciAir {},
    //     0
    // ));

    let pcs = Pcs::new(dft, val_mmcs, fri_config);

    type MyConfig = StarkConfig<Pcs, Challenge, Challenger>;
    let config = MyConfig::new(pcs);

    let mut challenger = Challenger::new(perm.clone());

    let proof = prove::<MyConfig, _>(&config, &air, &mut challenger, trace, &vec![]);

    std::fs::write(
        "proof_poseidon2_gl64.json",
        serde_json::to_string(&proof).unwrap(),
    )
    .unwrap();

    let mut challenger = Challenger::new(perm);
    verify(&config, &air, &mut challenger, &proof, &vec![]).unwrap();
}
