use ark_ff::{BigInteger, PrimeField};
use p3_field::{
    AbstractField,
    // PrimeField,
    PrimeField64,
};
use p3_goldilocks::{DiffusionMatrixGoldilocks, Goldilocks};
use p3_poseidon2::{
    DiffusionPermutation, MdsLightPermutation, Poseidon2, Poseidon2ExternalMatrixGeneral,
    M31_RC_16_30_U32,
};
use p3_symmetric::Permutation;
use rand::thread_rng;
use zkhash::fields::goldilocks::FpGoldiLocks;
use zkhash::poseidon::poseidon_instance_goldilocks::RC16;

fn goldilocks_from_ark_ff(input: FpGoldiLocks) -> Goldilocks {
    let as_bigint = input.into_bigint();
    let mut as_bytes = as_bigint.to_bytes_le();
    as_bytes.resize(8, 0);
    let as_u64 = u64::from_le_bytes(as_bytes[0..8].try_into().unwrap());
    Goldilocks::from_wrapped_u64(as_u64)
}

fn main() {
    const WIDTH: usize = 16;
    const D: u64 = 3;
    type Val = Goldilocks;

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

    let mut internal_constants = round_constants.clone()[4..26].iter().map(|x| x[0]).collect();

    let poseidon = Poseidon2::<Goldilocks, Poseidon2ExternalMatrixGeneral, DiffusionMatrixGoldilocks, WIDTH, D>::new(
        8,
        external_constants,
        external_linear_layer,
        22,
        internal_constants,
        internal_linear_layer
    );
    let input = [Goldilocks::zero(); WIDTH];
    let output = poseidon.permute(input);
    println!("output: {:?}", output);
}
