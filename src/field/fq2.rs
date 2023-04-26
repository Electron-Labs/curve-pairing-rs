use std::ops;
use num_bigint::BigInt;
use crate::field::fq::Fq;

#[derive(Debug, Clone)]
pub struct Fq2 {
    // non_residue refers to a non-quadratic residue in the base field Fq. A non-quadratic residue is an element of Fq that is not a square (i.e., not the square of any other element in Fq). It is used as a coefficient to compute the multiplication of elements in the field extension Fq2.
    pub non_residue: Fq,
    pub value: [Fq; 2],
    pub forbenius_coeff_c1: [Fq; 2]
}

impl Fq2 {
    pub fn fq2(value: [Fq; 2], non_residue: Fq, forbenius_coeff_c1: [Fq; 2]) -> Self {
        Self {
            non_residue,
            value: [value[0].clone(), value[1].clone()],
            forbenius_coeff_c1
        }
    }

    pub fn zero(non_residue: Fq, field_modulus: BigInt, forbenius_coeff_c1: [Fq; 2]) -> Self {
        return Self {
            non_residue,
            value: [Fq::zero(field_modulus.clone()), Fq::zero(field_modulus)],
            forbenius_coeff_c1
        }
    }

    pub fn one(non_residue: Fq, field_modulus: BigInt, forbenius_coeff_c1: [Fq; 2]) -> Self {
        return Self {
            non_residue,
            value: [Fq::one(field_modulus.clone()), Fq::zero(field_modulus)],
            forbenius_coeff_c1
        }
    }

    pub fn inverse(&self) -> Self {
         // High-Speed Software Implementation of the Optimal Ate Pairing over Barreto–Naehrig Curves .pdf
         // https://eprint.iacr.org/2010/354.pdf , algorithm 8
        let t0 = self.value[0].clone() * self.value[0].clone();
        let t1 = self.value[1].clone() * self.value[1].clone();
        let t2 = t0 - (t1 * self.non_residue.clone());
        let t3 = Fq::one(self.value[0].field_modulus.clone())/t2;
        let forbenius_coeff = self.forbenius_coeff_c1.clone();
        Self {
            non_residue: self.non_residue.clone(),
            value: [
                (self.value[0].clone() * t3.clone()),
                -(self.value[1].clone() * t3)
            ],
            forbenius_coeff_c1: forbenius_coeff
        }
    }

    pub fn forbenius_map(self, power: usize) -> Self {
        // Forbenius map is a field automorphism that maps each element of the field to its q-th power, where q is the order of the field.
        return Self {
            non_residue: self.non_residue,
            value: [
                self.value[0].clone(),
                self.forbenius_coeff_c1[power%2].clone() * self.value[1].clone()
            ],
            forbenius_coeff_c1: self.forbenius_coeff_c1
        }
    }

    pub fn mul_scalar(self, base: Fq) -> Self {
        let mut q = Fq2::zero(self.non_residue.clone(), self.value[0].clone().field_modulus, self.forbenius_coeff_c1.clone());
        let d = base.clone();
        let r = self.clone();
        let mut found_one = false;
        let d_bit_len = d.value.bits();
        for i in (0..d_bit_len).rev() {
            if found_one {
                q = q.clone() + q.clone();
            }
            if d.value.bit(i) {
                found_one = true;
                q = q+r.clone();
            }
        }
        return q;
    }
}

impl ops::Add<Fq2> for Fq2 {
    type Output = Fq2;
    fn add(self, rhs: Fq2) -> Self::Output {
        return Fq2::fq2(
            [
                    (self.value[0].clone()+rhs.value[0].clone()),
                    (self.value[1].clone()+rhs.value[1].clone())
                ],
            self.non_residue,
            self.forbenius_coeff_c1
        )
    }
}

impl ops::Sub<Fq2> for Fq2 {
    type Output = Fq2;
    fn sub(self, rhs: Fq2) -> Self::Output {
        return Fq2::fq2(
            [
                (self.value[0].clone()-rhs.value[0].clone()),
                (self.value[1].clone()-rhs.value[1].clone())
            ],
            self.non_residue,
            self.forbenius_coeff_c1
        )
    }
}

impl ops::Mul<Fq2> for Fq2 {
    // Multiplication and Squaring on Pairing-Friendly.pdf; Section 3 (Karatsuba)
    // https://pdfs.semanticscholar.org/3e01/de88d7428076b2547b60072088507d881bf1.pdf
    type Output = Fq2;
    fn mul(self, rhs: Fq2) -> Self::Output {
        let a0 = self.value[0].clone();
        let a1 = self.value[1].clone();
        let b0 = rhs.value[0].clone();
        let b1 = rhs.value[1].clone();

        let v0 = a0.clone() * b0.clone();
        let v1 = a1.clone() * b1.clone();

        let val1 = v0.clone() + (v1.clone() * self.non_residue.clone());
        let val2 = (a0.clone() + a1.clone()) * (b0.clone() + b1.clone()) - (v0.clone() + v1.clone());

        return Fq2::fq2(
            [val1, val2],
            self.non_residue,
            self.forbenius_coeff_c1
        )
    }
}

impl ops::Div<Fq2> for Fq2 {
    type Output = Fq2;
    fn div(self, rhs: Fq2) -> Self::Output {
        self * rhs.inverse()
    }
}

impl ops::Neg for Fq2 {
    type Output = Fq2;
    fn neg(self) -> Self::Output {
        Self {
            non_residue: self.non_residue,
            value: [-(self.value[0].clone()), -(self.value[1].clone())],
            forbenius_coeff_c1: self.forbenius_coeff_c1
        }
    }
}


impl PartialEq for Fq2 {
    fn eq(&self, other: &Self) -> bool {
        self.value[0] == other.value[0] && self.value[1] == other.value[1]
            && self.forbenius_coeff_c1[0] == other.forbenius_coeff_c1[0]
            && self.forbenius_coeff_c1[1] == other.forbenius_coeff_c1[1]
            && self.non_residue == other.non_residue
    }

    fn ne(&self, other: &Self) -> bool {
        !(self.value[0] == other.value[0] && self.value[1] == other.value[1]
            && self.forbenius_coeff_c1[0] == other.forbenius_coeff_c1[0]
            && self.forbenius_coeff_c1[1] == other.forbenius_coeff_c1[1]
            && self.non_residue == other.non_residue)
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;
    use super::*;

    #[test]
    fn test_fq2_zero() {
        let field_modulus: BigInt = BigInt::from(17);
        let non_residue = Fq::fq(field_modulus.clone(), BigInt::from(16));
        let forbenius_coeff_c1 = [Fq::fq(field_modulus.clone(),BigInt::from(1)),
            Fq::fq(field_modulus.clone(),BigInt::from_str("21888242871839275222246405745257275088696311157297823662689037894645226208582").unwrap())];
        let zero = Fq2::zero(non_residue.clone(), field_modulus.clone(), forbenius_coeff_c1.clone());
        assert_eq!(zero.value[0], Fq::zero(field_modulus.clone()));
        assert_eq!(zero.value[1], Fq::zero(field_modulus.clone()));
    }

    #[test]
    fn test_fq2_one() {
        let field_modulus: BigInt = BigInt::from(17);
        let non_residue = Fq::fq(field_modulus.clone(), BigInt::from(16));
        let forbenius_coeff_c1 = [Fq::fq(field_modulus.clone(),BigInt::from(1)),
            Fq::fq(field_modulus.clone(),BigInt::from_str("21888242871839275222246405745257275088696311157297823662689037894645226208582").unwrap())];
        let one = Fq2::one(non_residue, field_modulus.clone(), forbenius_coeff_c1);
        assert_eq!(one.value[0], Fq::one(field_modulus.clone()));
        assert_eq!(one.value[1], Fq::zero(field_modulus.clone()));
    }

    #[test]
    fn test_fq2_add() {
        let field_modulus = BigInt::from(7);
        let non_residue = Fq::fq(field_modulus.clone(), BigInt::from(-1));
        let fq_3 = Fq::fq(field_modulus.clone(), BigInt::from(3));
        let fq_4 = Fq::fq(field_modulus.clone(), BigInt::from(4));
        let forbenius_coeff_c1 = [Fq::fq(field_modulus.clone(),BigInt::from(1)),
            Fq::fq(field_modulus.clone(),BigInt::from_str("21888242871839275222246405745257275088696311157297823662689037894645226208582").unwrap())];
        let a = Fq2::fq2([fq_4.clone(), fq_4.clone()], non_residue.clone(), forbenius_coeff_c1.clone());
        let b = Fq2::fq2([fq_3.clone(), fq_4.clone()], non_residue.clone(), forbenius_coeff_c1.clone());
        let res = Fq2::fq2([Fq::zero(field_modulus.clone()), Fq::one(field_modulus.clone())], non_residue.clone(), forbenius_coeff_c1.clone());
        let c = a+b;
        assert_eq!(c, res);
    }

    #[test]
    fn test_fq2_sub() {
        let field_modulus = BigInt::from(7);
        let forbenius_coeff_c1 = [Fq::fq(field_modulus.clone(),BigInt::from(1)),
            Fq::fq(field_modulus.clone(),BigInt::from_str("21888242871839275222246405745257275088696311157297823662689037894645226208582").unwrap())];
        let non_residue = Fq::fq(field_modulus.clone(), BigInt::from(-1));
        let fq_5 = Fq::fq(field_modulus.clone(), BigInt::from(5));
        let fq_7 = Fq::fq(field_modulus.clone(), BigInt::from(7));
        let fq_3 = Fq::fq(field_modulus.clone(), BigInt::from(3));
        let fq_2 = Fq::fq(field_modulus.clone(), BigInt::from(2));
        let a = Fq2::fq2([fq_5.clone(), fq_3.clone()], non_residue.clone(), forbenius_coeff_c1.clone());
        let b = Fq2::fq2([fq_7.clone(), fq_2.clone()], non_residue.clone(), forbenius_coeff_c1.clone());
        let res = Fq2::fq2([fq_5.clone(), Fq::one(field_modulus.clone())],non_residue.clone(), forbenius_coeff_c1.clone());
        let c = a-b;
        assert_eq!(c, res);
    }

    #[test]
    fn test_fq2_mul() {
        let field_modulus = BigInt::from(7);
        let forbenius_coeff_c1 = [Fq::fq(field_modulus.clone(),BigInt::from(1)),
            Fq::fq(field_modulus.clone(),BigInt::from_str("21888242871839275222246405745257275088696311157297823662689037894645226208582").unwrap())];
        let non_residue = Fq::fq(field_modulus.clone(), BigInt::from(-1));
        let fq_4 = Fq::fq(field_modulus.clone(), BigInt::from(4));
        let fq_3 = Fq::fq(field_modulus.clone(), BigInt::from(3));
        let a = Fq2::fq2([fq_4.clone(), fq_4.clone()], non_residue.clone(), forbenius_coeff_c1.clone());
        let b = Fq2::fq2([fq_3.clone(), fq_4.clone()], non_residue.clone(), forbenius_coeff_c1.clone());
        let res = Fq2::fq2([fq_3.clone(), Fq::zero(field_modulus)], non_residue.clone(), forbenius_coeff_c1.clone());
        let c = a*b;
        assert_eq!(c,res);
    }

    #[test]
    fn test_fq2_neg() {
        let field_modulus = BigInt::from(17);
        let non_residue = Fq::fq(field_modulus.clone(), BigInt::from(16));
        let forbenius_coeff_c1 = [Fq::fq(field_modulus.clone(),BigInt::from(1)),
            Fq::fq(field_modulus.clone(),BigInt::from_str("21888242871839275222246405745257275088696311157297823662689037894645226208582").unwrap())];
        let fq_7 = Fq::fq(field_modulus.clone(), BigInt::from(7));
        let fq_8 = Fq::fq(field_modulus.clone(), BigInt::from(8));
        let fq_9 = Fq::fq(field_modulus.clone(), BigInt::from(9));
        let fq_10 = Fq::fq(field_modulus.clone(), BigInt::from(10));
        let a = Fq2::fq2([fq_9.clone(), fq_7.clone()], non_residue.clone(), forbenius_coeff_c1.clone());
        let a_neg = -a;
        let res = Fq2::fq2([fq_8.clone(), fq_10.clone()], non_residue.clone(), forbenius_coeff_c1.clone());
        assert_eq!(a_neg, res);
    }

    #[test]
    fn test_fq2_inverse() {
        let field_modulus = BigInt::from(7);
        let non_residue = Fq::fq(field_modulus.clone(), BigInt::from(-1));
        let forbenius_coeff_c1 = [Fq::fq(field_modulus.clone(),BigInt::from(1)),
            Fq::fq(field_modulus.clone(),BigInt::from_str("21888242871839275222246405745257275088696311157297823662689037894645226208582").unwrap())];
        let fq_4 = Fq::fq(field_modulus.clone(), BigInt::from(4));
        let fq_6 = Fq::fq(field_modulus.clone(), BigInt::from(6));
        let a = Fq2::fq2([fq_4.clone(), fq_4.clone()], non_residue.clone(), forbenius_coeff_c1.clone());
        let b = a.inverse();
        let d = a*b;
        println!("mult:{:?}", d);
        let res = Fq2::one(non_residue.clone(), field_modulus.clone(), forbenius_coeff_c1.clone());
        assert_eq!(d, res);
    }

    #[test]
    fn test_fq2_forbenius_map() {
        let field_modulus = BigInt::from_str("21888242871839275222246405745257275088696311157297823662689037894645226208583").unwrap();
        let non_residue = Fq::fq(field_modulus.clone(), BigInt::from_str("21888242871839275222246405745257275088696311157297823662689037894645226208582").unwrap());
        let forbenius_coeffs_c1 = [
            Fq::fq(field_modulus.clone(), BigInt::from(1)),
            Fq::fq(field_modulus.clone(), BigInt::from_str("21888242871839275222246405745257275088696311157297823662689037894645226208582").unwrap())
        ];
        let mut a = Fq2::fq2(
            [Fq::fq(field_modulus.clone(), BigInt::from(2)), Fq::fq(field_modulus.clone(), BigInt::from(3))],
            non_residue.clone(),
            forbenius_coeffs_c1.clone()
        );
        a = a.forbenius_map(3);
        let res = Fq2::fq2([Fq::fq(field_modulus.clone(), BigInt::from(2)),Fq::fq(field_modulus, BigInt::from_str("21888242871839275222246405745257275088696311157297823662689037894645226208580").unwrap())], non_residue, forbenius_coeffs_c1);
        assert_eq!(a, res)
    }

    #[test]
    fn test_fq2_mul_scalar() {
        let field_modulus = BigInt::from(17);
        let non_residue = Fq::fq(field_modulus.clone(), BigInt::from(16));
        let forbenius_coeffs_c1 = [Fq::fq(field_modulus.clone(),BigInt::from(1)),
            Fq::fq(field_modulus.clone(),BigInt::from_str("21888242871839275222246405745257275088696311157297823662689037894645226208582").unwrap())];
        let mut a = Fq2::fq2(
            [Fq::fq(field_modulus.clone(), BigInt::from(2)), Fq::fq(field_modulus.clone(), BigInt::from(3))],
            non_residue.clone(),
            forbenius_coeffs_c1.clone()
        );
        let b = Fq::fq(field_modulus.clone(), BigInt::from(3));
        let c = a.mul_scalar(b.clone());
        let res = Fq2::fq2([Fq::fq(field_modulus.clone(), BigInt::from(6)),Fq::fq(field_modulus, BigInt::from(9))], non_residue, forbenius_coeffs_c1);
        assert_eq!(c, res);
    }
}