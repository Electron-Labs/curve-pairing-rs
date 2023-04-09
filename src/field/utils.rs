use std::ops::{Div, Mul, Sub};
use ethnum::I256;
use num_bigint::BigInt;

fn egcd(a: BigInt, b: BigInt) -> (BigInt, BigInt, BigInt) {
    if a==BigInt::from(0) {
        return (b, BigInt::from(0), BigInt::from(1))
    } else {
        let (g, y, x) = egcd(b.clone()%a.clone(), a.clone());
        (g, x-(b.clone()/a.clone()) * y.clone(), y.clone())
    }
}

pub fn mod_inverse(a : BigInt, m: BigInt) -> BigInt {
    let a_ = a.clone();
    let m_ = m.clone();
    let (g, x, y) = egcd(a_, m_);
    assert_eq!(g, BigInt::from(1), "Modular inversion does not exist");
    (x % m)
}

