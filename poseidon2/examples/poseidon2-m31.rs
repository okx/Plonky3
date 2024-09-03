use p3_mersenne_31::{DiffusionMatrixMersenne31, Mersenne31};
use p3_field::{AbstractField, PrimeField, PrimeField64};
use p3_poseidon2::{
    DiffusionPermutation, MdsLightPermutation, Poseidon2, Poseidon2ExternalMatrixGeneral, RC_16_30_U32,
};
use p3_symmetric::Permutation;
use rand::thread_rng;


fn main() {
        const WIDTH: usize = 16;
    const D: u64 = 3;

    let RC_16_30_U32_M31: [[Mersenne31; WIDTH]; 30] = RC_16_30_U32
    .iter()
    .map(|round| round.map(Mersenne31::from_wrapped_u32))
    .collect::<Vec<_>>()
    .try_into()
    .unwrap();

    // Poseidon2ExternalMatrixGeneral, DiffusionMatrixMersenne31, 16, 5


    let mut rng = thread_rng();
    let external_linear_layer = Poseidon2ExternalMatrixGeneral::default();
    let internal_linear_layer = DiffusionMatrixMersenne31::default();

    let mut external_constants = RC_16_30_U32_M31.clone()[0..4].to_vec();
    external_constants.extend(RC_16_30_U32_M31.clone()[26..30].to_vec());
    // println!("external_constants: {:?}", external_constants);

    let mut internal_constants = RC_16_30_U32_M31.clone()[4..26].iter().map(|x| x[0]).collect();

    let poseidon = Poseidon2::<Mersenne31, Poseidon2ExternalMatrixGeneral, DiffusionMatrixMersenne31, WIDTH, D>::new(
        8,
        external_constants,
        external_linear_layer,
        22,
        internal_constants,
        internal_linear_layer
    );
    let input = [Mersenne31::zero(); WIDTH];
    let output = poseidon.permute(input);
    println!("output: {:?}", output);
    // let name = format!(
    //     "poseidon2::<{}, {}, {}>",
    //     type_name::<F::Packing>(),
    //     D,
    //     WIDTH
    // );
    // let id = BenchmarkId::new(name, WIDTH);
    // c.bench_with_input(id, &input, |b, &input| b.iter(|| ));
}