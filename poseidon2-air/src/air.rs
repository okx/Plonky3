///! the code is referenced from here: https://github.com/succinctlabs/sp1/pull/397
/// however, their implementation does not fit to plonky3; lots of parts have been tuned to make it fix plonky3
use alloc::vec::Vec;
use core::borrow::Borrow;

use p3_air::{Air, AirBuilder, BaseAir};
use p3_field::{AbstractField, Field};
use p3_matrix::Matrix;
use p3_mersenne_31::{POSEIDON2_INTERNAL_MATRIX_DIAG_16, POSEIDON2_INTERNAL_MATRIX_DIAG_16_SHIFTS};
use p3_poseidon2::{apply_mat4, M31_RC_16_30_U32};

use crate::columns::Poseidon2Cols;
use crate::{biguint_to_u64, num_cols, FullRound, PartialRound, SBox};

/// Assumes the field size is at least 16 bits.
///
/// ***WARNING***: this is a stub for now, not ready to use.
#[derive(Debug)]
pub struct Poseidon2Air<F: Field, const WIDTH: usize> {
    beginning_full_round_constants: [[F; 16]; 4],
    partial_round_constants: [F; 22],
    ending_full_round_constants: [[F; 16]; 4],
}

impl<F: Field, const WIDTH: usize> Poseidon2Air<F, WIDTH> {
    pub fn new() -> Self {
        let RC_16_30_U32_M31 = M31_RC_16_30_U32
            .iter()
            .map(|round| round.map(F::from_wrapped_u32))
            .collect::<Vec<_>>();
        // .try_into()
        // .unwrap();
        let beginning_full_round_constants: [[F; 16]; 4] =
            core::array::from_fn(|i| RC_16_30_U32_M31.get(i).unwrap().clone());
        let ending_full_round_constants: [[F; 16]; 4] =
            core::array::from_fn(|i| RC_16_30_U32_M31.get(i + 26).unwrap().clone());
        let partial_round_constants: [F; 22] =
            core::array::from_fn(|i| RC_16_30_U32_M31.get(i + 4).unwrap().clone()[0]);

        Self {
            beginning_full_round_constants: beginning_full_round_constants,
            partial_round_constants: partial_round_constants,
            ending_full_round_constants: ending_full_round_constants,
        }
    }
}

impl<F: Field, const WIDTH: usize> BaseAir<F> for Poseidon2Air<F, WIDTH> {
    fn width(&self) -> usize {
        num_cols::<WIDTH>()
    }
}

impl<AB: AirBuilder, const WIDTH: usize> Air<AB> for Poseidon2Air<AB::F, WIDTH> {
    #[inline]
    fn eval(&self, builder: &mut AB) {
        let main = builder.main();
        let local = main.row_slice(0);
        let local: &Poseidon2Cols<AB::Var, WIDTH> = (*local).borrow();

        //TODO: this is fixed; but it might depends on the WIDTH and degree of S BOX
        let rounds_f = 8;
        let rounds_p = 22;
        let rounds = rounds_f + rounds_p;

        // Convert the u32 round constants to field elements.
        // [[AB::F; WIDTH]; 30]
        let constants = M31_RC_16_30_U32
            .iter()
            .map(|round| round.map(AB::F::from_wrapped_u32))
            .collect::<Vec<_>>();

        // Apply the round constants.
        //
        // Initial Layer: Don't apply the round constants.
        // External Layers: Apply the round constants.
        // Internal Layers: Only apply the round constants to the first element.
        for i in 0..WIDTH {
            let mut result: AB::Expr = local.input[i].into();
            for r in 0..rounds {
                if i == 0 {
                    result += local.rounds[r + 1]
                        * constants[r][i]
                        * (local.is_external + local.is_internal);
                } else {
                    result += local.rounds[r + 1] * constants[r][i] * local.is_external;
                }
            }
            builder.assert_eq(result, local.add_rc[i]);
        }

        // Apply the sbox.
        //
        // To differentiate between external and internal layers, we use a masking operation
        // to only apply the state change to the first element for internal layers.
        for i in 0..WIDTH {
            let sbox_deg_3 = local.add_rc[i] * local.add_rc[i] * local.add_rc[i];
            builder.assert_eq(sbox_deg_3, local.sbox_deg_3[i]);
            // let sbox_deg_7 = local.sbox_deg_3[i] * local.sbox_deg_3[i] * local.add_rc[i];
            // builder.assert_eq(sbox_deg_7, local.sbox_deg_7[i]);
        }
        let sbox_result: [AB::Expr; WIDTH] = local
            .sbox_deg_3
            .iter()
            .enumerate()
            .map(|(i, x)| {
                // The masked first result of the sbox.
                //
                // Initial Layer: Pass through the result of the round constant layer.
                // External Layer: Pass through the result of the sbox layer.
                // Internal Layer: Pass through the result of the sbox layer.
                if i == 0 {
                    local.is_initial * local.add_rc[i] + (AB::Expr::one() - local.is_initial) * *x
                }
                // The masked result of the rest of the sbox.
                //
                // Initial layer: Pass through the result of the round constant layer.
                // External layer: Pass through the result of the sbox layer.
                // Internal layer: Pass through the result of the round constant layer.
                else {
                    (local.is_initial + local.is_internal) * local.add_rc[i]
                        + (AB::Expr::one() - (local.is_initial + local.is_internal)) * *x
                }
            })
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();

        // EXTERNAL LAYER + INITIAL LAYER
        {
            // First, we apply M_4 to each consecutive four elements of the state.
            // In Appendix B's terminology, this replaces each x_i with x_i'.
            let mut state: [AB::Expr; WIDTH] = sbox_result.clone();
            for i in (0..WIDTH).step_by(4) {
                let mut mat4_state = core::array::from_fn(|k| state[i + k].clone());
                apply_mat4(&mut mat4_state);
                state[i..i + 4].clone_from_slice(&mat4_state);
            }

            // Now, we apply the outer circulant matrix (to compute the y_i values).
            //
            // We first precompute the four sums of every four elements.
            let sums: [AB::Expr; 4] = core::array::from_fn(|k| {
                (0..WIDTH)
                    .step_by(4)
                    .map(|j| state[j + k].clone())
                    .sum::<AB::Expr>()
            });

            // The formula for each y_i involves 2x_i' term and x_j' terms for each j that equals i mod 4.
            // In other words, we can add a single copy of x_i' to the appropriate one of our precomputed sums.
            for i in 0..WIDTH {
                state[i] += sums[i % 4].clone();
                builder
                    .when(local.is_external + local.is_initial)
                    .assert_eq(state[i].clone(), local.output[i]);
            }
        }

        // INTERNAL LAYER
        {
            // Use a simple matrix multiplication as the permutation.
            // let mut state: [AB::Expr; WIDTH] = sbox_result.clone();
            // let matmul_constants: [<<AB as AirBuilder>::Expr as AbstractField>::F; WIDTH] =
            //     MATRIX_DIAG_16_M31_U32
            //         .iter()
            //         .map(|x| <<AB as AirBuilder>::Expr as AbstractField>::F::from_wrapped_u32(*x))
            //         .collect::<Vec<_>>()
            //         .try_into()
            //         .unwrap();
            // matmul_internal(&mut state, matmul_constants);
            let mut state: [AB::Expr; WIDTH] = sbox_result.clone();
            let part_sum: AB::Expr = state.clone()
                .iter()
                .skip(1)
                .map(|x| {
                    x.clone()
                })
                .sum();
            let full_sum = part_sum.clone() + state[0].clone();
            let s0 = part_sum.clone() + (-state[0].clone());
            // state[0] = F::from_canonical_u32(biguint_to_u32(from_u62(s0).as_canonical_biguint()));
            state[0] = s0;
            for i in 1..16 {
                let si = full_sum.clone()
                    + (state[i].clone()
                        * (AB::F::from_canonical_u32(
                            1 << POSEIDON2_INTERNAL_MATRIX_DIAG_16_SHIFTS[i-1],
                        )));
                state[i] = si;
            }

            for i in 0..WIDTH {
                builder
                    .when(local.is_internal)
                    .assert_eq(state[i].clone(), local.output[i]);
            }
        }

        // Range check all flags.
        for i in 0..local.rounds.len() {
            builder.assert_bool(local.rounds[i]);
        }
        builder.assert_bool(local.is_initial);
        builder.assert_bool(local.is_external);
        builder.assert_bool(local.is_internal);
        builder.assert_bool(local.is_initial + local.is_external + local.is_internal);

        // Constrain the initial flag.
        builder.assert_eq(local.is_initial, local.rounds[0]);

        // Constrain the external flag.
        let is_external_first_half = (0..4).map(|i| local.rounds[i + 1].into()).sum::<AB::Expr>();
        let is_external_second_half = (26..30)
            .map(|i| local.rounds[i + 1].into())
            .sum::<AB::Expr>();
        builder.assert_eq(
            local.is_external,
            is_external_first_half + is_external_second_half,
        );

        // Constrain the internal flag.
        let is_internal = (4..26)
            .map(|i| local.rounds[i + 1].into())
            .sum::<AB::Expr>();
        builder.assert_eq(local.is_internal, is_internal);
    }
}


#[inline]
fn eval_full_round<
    AB: AirBuilder,
    const WIDTH: usize,
    const SBOX_DEGREE: usize,
    const SBOX_REGISTERS: usize,
>(
    state: &mut [AB::Expr; WIDTH],
    full_round: &FullRound<AB::Var, WIDTH, SBOX_DEGREE, SBOX_REGISTERS>,
    round_constants: &[AB::F; WIDTH],
    builder: &mut AB,
) {
    for (i, (s, r)) in state.iter_mut().zip(round_constants.iter()).enumerate() {
        *s = s.clone() + *r;
        eval_sbox(&full_round.sbox[i], s, builder);
    }
    // L::matmul_external(state);
}

#[inline]
fn eval_partial_round<
    AB: AirBuilder,
    const WIDTH: usize,
    const SBOX_DEGREE: usize,
    const SBOX_REGISTERS: usize,
>(
    state: &mut [AB::Expr; WIDTH],
    partial_round: &PartialRound<AB::Var, WIDTH, SBOX_DEGREE, SBOX_REGISTERS>,
    round_constant: &AB::F,
    builder: &mut AB,
) {
    state[0] = state[0].clone() + *round_constant;
    eval_sbox(&partial_round.sbox, &mut state[0], builder);
    // L::matmul_internal(state, internal_matrix_diagonal);
}

/// Evaluates the S-BOX over a degree-`1` expression `x`.
///
/// # Panics
///
/// This method panics if the number of `REGISTERS` is not chosen optimally for the given
/// `DEGREE` or if the `DEGREE` is not supported by the S-BOX. The supported degrees are
/// `3`, `5`, `7`, and `11`.
///
/// # Efficiency Note
///
/// This method computes the S-BOX by computing the cube of `x` and then successively
/// multiplying the running sum by the cube of `x` until the last multiplication where we use
/// the appropriate power to reach the final product:
///
/// ```text
/// (x^3) * (x^3) * ... * (x^k) where k = d mod 3
/// ```
///
/// The intermediate powers are stored in the auxiliary column registers. To maximize the
/// efficiency of the registers we try to do three multiplications per round. This algorithm
/// only multiplies the cube of `x` but a more optimal product would be to find the base-3
/// decomposition of the `DEGREE` and use that to generate the addition chain. Even this is not
/// the optimal number of multiplications for all possible degrees, but for the S-BOX powers we
/// are interested in for Poseidon2 (namely `3`, `5`, `7`, and `11`), we get the optimal number
/// with this algorithm. We use the following register table:
///
/// | `DEGREE` | `REGISTERS` |
/// |:--------:|:-----------:|
/// | `3`      | `1`         |
/// | `5`      | `2`         |
/// | `7`      | `3`         |
/// | `11`     | `3`         |
///
/// We record this table in [`Self::OPTIMAL_REGISTER_COUNT`] and this choice of registers is
/// enforced by this method.
#[inline]
fn eval_sbox<AB, const DEGREE: usize, const REGISTERS: usize>(
    sbox: &SBox<AB::Var, DEGREE, REGISTERS>,
    x: &mut AB::Expr,
    builder: &mut AB,
) where
    AB: AirBuilder,
{
    // assert_ne!(REGISTERS, 0, "The number of REGISTERS must be positive.");
    // assert!(DEGREE <= 11, "The DEGREE must be less than or equal to 11.");
    // assert_eq!(
    //     REGISTERS,
    //     Self::OPTIMAL_REGISTER_COUNT[DEGREE],
    //     "The number of REGISTERS must be optimal for the given DEGREE."
    // );

    let x2 = x.square();
    let x3 = x2.clone() * x.clone();

    load(sbox, 0, x3.clone(), builder);
    if REGISTERS == 1 {
        // println!("x3: {:?}", x3);
        *x = sbox.0[0].into();
        // println!("sbox: {:?}", sbox.0[0].into());
        return;
    }
    if DEGREE == 11 {
        (1..REGISTERS - 1).for_each(|j| load_product(sbox, j, &[0, 0, j - 1], builder));
    } else {
        (1..REGISTERS - 1).for_each(|j| load_product(sbox, j, &[0, j - 1], builder));
    }
    load_last_product(sbox, x.clone(), x2, x3, builder);
    *x = sbox.0[REGISTERS - 1].into();
}

/// Loads `value` into the `i`-th S-BOX register.
#[inline]
fn load<AB, const SBOX_DEGREE: usize, const SBOX_REGISTERS: usize>(
    _sbox: &SBox<AB::Var, SBOX_DEGREE, SBOX_REGISTERS>,
    _i: usize,
    _value: AB::Expr,
    _builder: &mut AB,
) where
    AB: AirBuilder,
{
    // println!("load sbox i: {:}, value: {:?}", _i, _value);
    _builder.assert_eq(_sbox.0[_i].into(), _value);
}

/// Loads the product over all `product` indices the into the `i`-th S-BOX register.
#[inline]
fn load_product<AB, const SBOX_DEGREE: usize, const SBOX_REGISTERS: usize>(
    _sbox: &SBox<AB::Var, SBOX_DEGREE, SBOX_REGISTERS>,
    _i: usize,
    _product: &[usize],
    _builder: &mut AB,
) where
    AB: AirBuilder,
{
    // assert!(
    //     product.len() <= 3,
    //     "Product is too big. We can only compute at most degree-3 constraints."
    // );
    // load(
    //     sbox,
    //     i,
    //     product.iter().map(|j| AB::Expr::from(self.0[*j])).product(),
    //     builder,
    // );
}

/// Loads the final product into the last S-BOX register. The final term in the product is
/// `pow(x, DEGREE % 3)`.
#[inline]
fn load_last_product<AB, const SBOX_DEGREE: usize, const SBOX_REGISTERS: usize>(
    _sbox: &SBox<AB::Var, SBOX_DEGREE, SBOX_REGISTERS>,
    _x: AB::Expr,
    _x2: AB::Expr,
    _x3: AB::Expr,
    _builder: &mut AB,
) where
    AB: AirBuilder,
{
    // load(
    //     sbox,
    //     REGISTERS - 1,
    //     [x3, x, x2][DEGREE % 3].clone() * AB::Expr::from(self.0[REGISTERS - 2]),
    //     builder,
    // );
}
