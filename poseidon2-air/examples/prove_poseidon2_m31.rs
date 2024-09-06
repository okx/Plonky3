use std::fmt::Debug;
use std::marker::PhantomData;

use p3_challenger::{DuplexChallenger, HashChallenger, SerializingChallenger32};
use p3_circle::CirclePcs;
use p3_commit::ExtensionMmcs;
use p3_dft::Radix2DitParallel;
use p3_field::extension::BinomialExtensionField;
use p3_field::{AbstractField, Field};
use p3_fri::{FriConfig, TwoAdicFriPcs};
use p3_keccak::Keccak256Hash;
use p3_merkle_tree::FieldMerkleTreeMmcs;
use p3_mersenne_31::{DiffusionMatrixMersenne31, Mersenne31};
use p3_poseidon2::{Poseidon2, Poseidon2ExternalMatrixGeneral, M31_RC_16_30_U32};
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

const WIDTH: usize = 16;
const SBOX_DEGREE: usize = 3;
const SBOX_REGISTERS: usize = 1;
const HALF_FULL_ROUNDS: usize = 4;
const PARTIAL_ROUNDS: usize = 20;

const NUM_HASHES: usize = 1 << 4;

// F::from_wrapped_u32(M31_RC_16_30_U32[r - 1][j]

fn main() {
    let env_filter = EnvFilter::builder()
        .with_default_directive(LevelFilter::INFO.into())
        .from_env_lossy();

    Registry::default()
        .with(env_filter)
        .with(ForestLayer::default())
        .init();

    type Val = Mersenne31;
    type Challenge = Val;
    type ByteHash = Keccak256Hash;
    type FieldHash = SerializingHasher32<ByteHash>;
    let byte_hash = ByteHash {};
    let field_hash = FieldHash::new(Keccak256Hash {});

    type MyCompress = CompressionFunctionFromHasher<u8, ByteHash, 2, 32>;
    let compress = MyCompress::new(byte_hash);

    type ValMmcs = FieldMerkleTreeMmcs<Val, u8, FieldHash, MyCompress, 32>;
    let val_mmcs = ValMmcs::new(field_hash, compress);

    type ChallengeMmcs = ExtensionMmcs<Val, Challenge, ValMmcs>;
    let challenge_mmcs = ChallengeMmcs::new(val_mmcs.clone());

    type Dft = Radix2DitParallel;
    let dft = Dft {};

    type Challenger = SerializingChallenger32<Val, HashChallenger<u8, ByteHash, 32>>;

    let external_linear_layer = Poseidon2ExternalMatrixGeneral::default();
    let internal_linear_layer = DiffusionMatrixMersenne31::default();
    let air: Poseidon2Air<Val, WIDTH> = Poseidon2Air::new();
    // let mut inputs = (0..NUM_HASHES).map(|i| core::array::from_fn(|j| Mersenne31::from_canonical_u32(0))).collect::<Vec<_>>();
    let mut input = core::array::from_fn(|j| Mersenne31::from_canonical_u32(0));

    let rc: Vec<[Val; WIDTH]> = M31_RC_16_30_U32
        .iter()
        .map(|x| core::array::from_fn(|i| Val::from_wrapped_u32(x[i])))
        .collect();
    let trace = generate_trace::<Val, WIDTH>(&mut input, rc, FieldType::M31);

    let fri_config = FriConfig {
        log_blowup: 1,
        num_queries: 100,
        proof_of_work_bits: 16,
        mmcs: challenge_mmcs,
    };

    type Pcs = CirclePcs<Val, ValMmcs, ChallengeMmcs>;
    let pcs = Pcs {
        mmcs: val_mmcs,
        fri_config,
        _phantom: PhantomData,
    };
    type MyConfig = StarkConfig<Pcs, Challenge, Challenger>;
    let config = MyConfig::new(pcs);

    let mut challenger = Challenger::from_hasher(vec![], byte_hash);
    let proof = prove(&config, &air, &mut challenger, trace, &vec![]);

    std::fs::write(
        "proof_poseidon2_m31.json",
        serde_json::to_string(&proof).unwrap(),
    )
    .unwrap();

    // let mut challenger = Challenger::new(perm);
    verify(&config, &air, &mut challenger, &proof, &vec![]);

    //
}
