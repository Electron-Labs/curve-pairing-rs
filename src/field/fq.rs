use std::ops;
use crate::field::utils::mod_inverse;

#[derive(Debug, Clone)]
struct Fq {
    pub field_modulus: i128,
    pub value: i128
}

impl Fq {
    pub fn fq(field_modulus: i128, value: i128) -> Self {
        let val = (value % field_modulus);
        Self { field_modulus, value: val }
    }

    pub fn zero(field_modulus: i128) -> Self {
        Self { field_modulus, value: 0 }
    }

    pub fn one(field_modulus: i128) -> Self {
        Self { field_modulus, value: 1}
    }
}

impl ops::Add<Fq> for Fq {
    type Output = Fq;
    fn add(self, rhs: Fq) -> Self::Output {
        assert_eq!(self.field_modulus, rhs.field_modulus, "Field modulus of both Fq mismatch");
        println!("Fq - ({:?}) + Fq - ({:?})", self, rhs);
        Fq::fq(self.field_modulus, (self.value + rhs.value) % self.field_modulus)
    }
}

impl ops::Sub<Fq> for Fq {
    type Output = Fq;
    fn sub(self, rhs: Fq) -> Self::Output {
        assert_eq!(self.field_modulus, rhs.field_modulus, "Field modulus of both Fq mismatch");
        println!("Fq - ({:?}) - Fq - ({:?})", self, rhs);
        Fq::fq(self.field_modulus, (self.value - rhs.value) % self.field_modulus)
    }
}

impl ops::Mul<Fq> for Fq {
    type Output = Fq;
    fn mul(self, rhs: Fq) -> Self::Output {
        assert_eq!(self.field_modulus, rhs.field_modulus, "Field modulus of both Fq mismatch");
        println!("Fq - ({:?}) * Fq - ({:?})", self, rhs);
        Fq::fq(self.field_modulus, (self.value * rhs.value) % self.field_modulus)
    }
}

impl ops::Div<Fq> for Fq {
    type Output = Fq;
    fn div(self, rhs: Fq) -> Self::Output {
        assert_eq!(self.field_modulus, rhs.field_modulus, "Field modulus of both Fq mismatch");
        println!("Fq - ({:?}) / Fq - ({:?})", self, rhs);
        Fq::fq(self.field_modulus, (self.value * mod_inverse(rhs.value, self.field_modulus))%self.field_modulus)
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
    const field_modulus: i128 = 100_000_000;

    #[test]
    fn test_fq_zero() {
        let zero_auto = Fq::zero(field_modulus);
        let zero_custom = Fq::fq(field_modulus, 0);
        assert_eq!(zero_auto, zero_custom);
    }

    #[test]
    fn test_fq_one() {
        let one_auto = Fq::one(field_modulus);
        let one_custom = Fq::fq(field_modulus, 1);
        assert_eq!(one_auto, one_custom);
    }

    #[test]
    fn test_fq_add() {
        let a = Fq::zero(field_modulus);
        let b = Fq::one(field_modulus);
        let c = a+ b.clone();
        assert_eq!(c, b);
    }

    #[test]
    fn test_fq_sub() {
        let a = Fq::one(field_modulus);
        let b = Fq::one(field_modulus);
        let c = Fq::zero(field_modulus);
        assert_eq!(a-b, c);
    }

    #[test]
    fn test_fq_mul() {
        let a = Fq::fq(field_modulus, 3);
        let b = Fq::fq(field_modulus, 2);
        let c = Fq::fq(field_modulus, 6);
        assert_eq!((a*b), c);
    }

    #[test]
    fn test_fq_div() {
        let fm = 17;
        let a = Fq::fq(fm, 2);
        let b = Fq::fq(fm, 9);
        let res = Fq::fq(fm ,4);
        assert_eq!(res, a/b);
    }
}

