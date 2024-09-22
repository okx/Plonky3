//! And AIR for the Poseidon2 permutation.

// #![no_std]

extern crate alloc;

mod air;
mod columns;
mod generation;

pub use air::*;
pub use columns::*;
pub use generation::*;

use p3_goldilocks::{DiffusionMatrixGoldilocks, Goldilocks};
use zkhash::fields::goldilocks::FpGoldiLocks;
use ark_ff::{BigInteger, PrimeField};
use p3_field::{AbstractField, Field};

pub fn goldilocks_from_ark_ff(input: FpGoldiLocks) -> Goldilocks {
    let as_bigint = input.into_bigint();
    let mut as_bytes = as_bigint.to_bytes_le();
    as_bytes.resize(8, 0);
    let as_u64 = u64::from_le_bytes(as_bytes[0..8].try_into().unwrap());
    Goldilocks::from_wrapped_u64(as_u64)
}