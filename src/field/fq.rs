use std::ops;
use crate::field::utils::mod_inverse;
use ethnum::I256;
use num_bigint::BigInt;

#[derive(Debug, Clone)]
pub struct Fq {
    pub field_modulus: BigInt,
    pub value: BigInt
}

impl Fq{
    pub fn fq(field_modulus: BigInt, value: BigInt) -> Self {
        let fm_ = field_modulus.clone();
        if value< BigInt::from(0) {
            Self { field_modulus: fm_.clone(), value: (value + fm_.clone())%(fm_.clone()) }
        } else {
            Self { field_modulus: fm_.clone(), value: value%(fm_.clone()) }
        }
    }

    pub fn zero(field_modulus: BigInt) -> Self {
        Self { field_modulus, value: BigInt::from(0) }
    }

    pub fn one(field_modulus: BigInt) -> Self {
        Self { field_modulus, value: BigInt::from(1) }
    }
}

impl ops::Add<Fq> for Fq {
    type Output = Fq;
    fn add(self, rhs: Fq) -> Self::Output {
        assert_eq!(self.field_modulus, rhs.field_modulus, "Field modulus of both Fq mismatch");
        // println!("Fq - ({:?}) + Fq - ({:?})", self, rhs);
        Fq::fq(self.field_modulus.clone(), (self.value + rhs.value) % self.field_modulus.clone())
    }
}

impl ops::Sub<Fq> for Fq {
    type Output = Fq;
    fn sub(self, rhs: Fq) -> Self::Output {
        assert_eq!(self.field_modulus, rhs.field_modulus, "Field modulus of both Fq mismatch");
        // println!("Fq - ({:?}) - Fq - ({:?})", self, rhs);
        Fq::fq(self.field_modulus.clone(), (self.value - rhs.value) % self.field_modulus.clone())
    }
}

impl ops::Mul<Fq> for Fq {
    type Output = Fq;
    fn mul(self, rhs: Fq) -> Self::Output {
        assert_eq!(self.field_modulus, rhs.field_modulus, "Field modulus of both Fq mismatch");
        // println!("Fq - ({:?}) * Fq - ({:?})", self, rhs);
        Fq::fq(self.field_modulus.clone(), (self.value * rhs.value) % self.field_modulus.clone())
    }
}

impl ops::Div<Fq> for Fq {
    type Output = Fq;
    fn div(self, rhs: Fq) -> Self::Output {
        assert_eq!(self.field_modulus, rhs.field_modulus, "Field modulus of both Fq mismatch");
        // println!("Fq - ({:?}) / Fq - ({:?})", self, rhs);
        Fq::fq(self.field_modulus.clone(), (self.value * mod_inverse(rhs.value.clone(), self.field_modulus.clone()))%self.field_modulus.clone())
    }
}

impl ops::Neg for Fq {
    type Output = Fq;
    fn neg(self) -> Self::Output {
        // Self { field_modulus: self.field_modulus, value: I256::from(-1) * self.value }
        Fq::fq(self.field_modulus, (BigInt::from(-1) * self.value))
    }
}

impl PartialEq for Fq {
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value && self.field_modulus == other.field_modulus
    }

    fn ne(&self, other: &Self) -> bool {
        self.value != other.value || self.field_modulus != other.field_modulus
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const field_modulus: i128 = 100000000;

    #[test]
    fn test_fq_zero() {
        let zero_auto = Fq::zero(BigInt::from(field_modulus));
        let zero_custom = Fq::fq(BigInt::from(field_modulus), BigInt::from(0));
        assert_eq!(zero_auto, zero_custom);
    }

    #[test]
    fn test_fq_one() {
        let one_auto = Fq::one(BigInt::from(field_modulus));
        let one_custom = Fq::fq(BigInt::from(field_modulus), BigInt::from(1));
        assert_eq!(one_auto, one_custom);
    }

    #[test]
    fn test_fq_add() {
        let a = Fq::zero(BigInt::from(field_modulus));
        let b = Fq::one(BigInt::from(field_modulus));
        let c = a+ b.clone();
        assert_eq!(c, b);
    }

    #[test]
    fn test_fq_sub() {
        let a = Fq::one(BigInt::from(field_modulus));
        let b = Fq::one(BigInt::from(field_modulus));
        let c = Fq::zero(BigInt::from(field_modulus));
        assert_eq!(a-b, c);
    }

    #[test]
    fn test_fq_mul() {
        let a = Fq::fq(BigInt::from(field_modulus), BigInt::from(3));
        let b = Fq::fq(BigInt::from(field_modulus), BigInt::from(2));
        let c = Fq::fq(BigInt::from(field_modulus), BigInt::from(6));
        assert_eq!((a*b), c);
    }

    #[test]
    fn test_fq_div() {
        let fm = BigInt::from(17);
        let a = Fq::fq(fm.clone(), BigInt::from(2));
        let b = Fq::fq(fm.clone(), BigInt::from(9));
        let res = Fq::fq(fm ,BigInt::from(4));
        assert_eq!(res, a/b);
    }
}

