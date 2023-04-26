use std::ops;
use num_bigint::BigInt;
use crate::field::fq::Fq;
use crate::field::fq2::Fq2;


#[derive(Debug, Clone)]
pub struct G2 {
    pub x: Fq2,
    pub y: Fq2,
    pub z: Fq2
}

impl G2 {
    pub fn new(x: Fq2, y:Fq2, z:Fq2) -> Self {
        Self {
            x, y, z
        }
    }

    pub fn isZero(&self) -> bool {
        return self.z == Fq2::zero(self.z.clone().non_residue, self.z.value[0].clone().field_modulus, self.z.clone().forbenius_coeff_c1)
    }

    pub fn double(&self) -> Self {
        if self.isZero() {
            return self.clone()
        }
        let x = self.x.clone();
        let y = self.y.clone();
        let z = self.z.clone();

        let a = x.clone() * x.clone();
        let b = y.clone() * y.clone();
        let c = b.clone() * b.clone();

        let t0 = x.clone() + b.clone();
        let t1 = t0.clone() * t0.clone();
        let t2 = t1.clone() - a.clone();
        let t3 = t2.clone() - c.clone();

        let d = t3.clone() + t3.clone();
        let e = a.clone() + a.clone() + a.clone();
        let f = e.clone() * e.clone();

        let t4 = d.clone() + d.clone();
        let x3 = f.clone() - t4.clone();

        let t5 = d.clone() - x3.clone();
        let c_2 = c.clone() + c.clone();
        let c_4 = c_2.clone() + c_2.clone();
        let t6 = c_4.clone() + c_4.clone();
        let t7 = e.clone() * t5.clone();
        let y3 = t7.clone() - t6.clone();

        let t8 = y.clone() * z.clone();
        let z3 = t8.clone() + t8.clone();

        return G2::new(x3, y3, z3)
    }

    pub fn affine(&self) -> Self {
        let non_residue = self.x.clone().non_residue;
        let field_modulus = self.x.value[0].clone().field_modulus;
        let forbenius_coeffs = self.x.clone().forbenius_coeff_c1;
        if self.isZero() {
            return G2::new(Fq2::zero(non_residue.clone(), field_modulus.clone(), forbenius_coeffs.clone()),
                           Fq2::zero(non_residue.clone(), field_modulus.clone(), forbenius_coeffs.clone()),
                           Fq2::zero(non_residue.clone(), field_modulus.clone(), forbenius_coeffs.clone()))
        }
        let zinv = Fq2::one(non_residue.clone(), field_modulus.clone(), forbenius_coeffs.clone()) / self.z.clone();
        let zinv2 = zinv.clone() * zinv.clone();
        let x = self.x.clone() * zinv2.clone();

        let zinv3 = zinv2.clone() * zinv.clone();
        let y = self.y.clone() * zinv3.clone();
        return G2::new(x, y, Fq2::one(non_residue.clone(), field_modulus.clone(), forbenius_coeffs.clone()))
    }

    pub fn mul_scalar(&self, base: Fq) -> Self {
        let non_residue = self.x.clone().non_residue;
        let field_modulus = self.x.value[0].clone().field_modulus;
        let forbenius_coeffs = self.x.clone().forbenius_coeff_c1;
        let mut q = G2::new(Fq2::zero(non_residue.clone(), field_modulus.clone(), forbenius_coeffs.clone()),
                        Fq2::zero(non_residue.clone(), field_modulus.clone(), forbenius_coeffs.clone()),
                        Fq2::zero(non_residue.clone(), field_modulus.clone(), forbenius_coeffs.clone()));
        let d = base.clone();
        let r = self.clone();
        let mut found_one = false;

        let d_bit_len = d.value.bits();
        for i in (0..d_bit_len).rev() {
            if found_one {
                q = q.double();
            }
            if d.value.bit(i) {
                found_one = true;
                q = q+r.clone();
            }
        }
        return q;
    }
}

impl ops::Add<G2> for G2 {
    type Output = G2;

    fn add(self, rhs: G2) -> Self::Output {
        if self.isZero() {
            return rhs.clone()
        }
        if rhs.isZero() {
            return self.clone()
        }
        let x1 = self.x.clone();
        let y1 = self.y.clone();
        let z1 = self.z.clone();
        let x2 = rhs.x.clone();
        let y2 = rhs.y.clone();
        let z2 = rhs.z.clone();

        let z1_2 = z1.clone()*z1.clone();
        let z2_2 = z2.clone()*z2.clone();

        let u1 = x1.clone() * z2_2.clone();
        let u2 = x2.clone() * z1_2.clone();

        let t0 = z2.clone() * z2_2.clone();
        let s1 = y1.clone() * t0.clone();

        let t1 = z1.clone() * z1_2.clone();
        let s2 = y2.clone() * t1.clone();

        let h = u2.clone() - u1.clone();
        let t2 = h.clone() + h.clone();
        let i = t2.clone() * t2.clone();
        let j = h.clone() * i.clone();
        let t3 = s2.clone() - s1.clone();
        let r = t3.clone() + t3.clone();
        let v = u1.clone() * i.clone();
        let t4 = r.clone() * r.clone();
        let t5 = v.clone() + v.clone();
        let t6 = t4.clone() - j.clone();
        let x3 = t6.clone() - t5.clone();
        let t7 = v.clone() - x3.clone();
        let t8 = s1.clone() * j.clone();
        let t9 = t8.clone() + t8.clone();
        let t10 = r.clone() * t7.clone();

        let y3 = t10.clone() - t9.clone();

        let t11 = z1.clone() + z2.clone();
        let t12 = t11.clone() * t11.clone();
        let t13 = t12.clone() - z1_2.clone();
        let t14 = t13.clone() - z2_2.clone();
        let z3 = t14.clone() * h.clone();

        return G2::new(x3, y3, z3)
    }
}

impl ops::Neg for G2 {
    type Output = G2;

    fn neg(self) -> Self::Output {
        return G2::new(self.x, -self.y, self.z)
    }
}

impl ops::Sub for G2 {
    type Output = G2;

    fn sub(self, rhs: Self) -> Self::Output {
        return self + (-rhs)
    }
}

impl PartialEq for G2 {
    fn eq(&self, other: &Self) -> bool {
        if self.isZero() {
            return other.isZero()
        }
        if other.isZero() {
            return self.isZero()
        }

        let z1_2 = self.z.clone() * self.z.clone();
        let z2_2 = other.z.clone() * other.z.clone();

        let u1 = self.x.clone() * z2_2.clone();
        let u2 = other.x.clone() * z1_2.clone();

        let z1_3 = self.z.clone() * z1_2.clone();
        let z2_3 = other.z.clone() * z2_2.clone();

        let s1 = self.y.clone() * z2_3.clone();
        let s2 = other.y.clone() * z1_3.clone();

        return u1 == u2 && s1 == s2;
    }
}

#[cfg(test)]
mod tests{
    use std::str::FromStr;
    use super::*;

    fn get_fq(val: &str) -> Fq {
        let val_ = BigInt::from_str(val).unwrap();
        let field_modulus = BigInt::from_str("21888242871839275222246405745257275088696311157297823662689037894645226208583").unwrap();
        return Fq::fq(field_modulus, val_);
    }

    #[test]
    fn test_g2_add() {
        let field_modulus = BigInt::from_str("21888242871839275222246405745257275088696311157297823662689037894645226208583").unwrap();
        let non_residue = Fq::fq(field_modulus.clone(), BigInt::from_str("21888242871839275222246405745257275088696311157297823662689037894645226208582").unwrap());
        let forbenius_coeffs_c1 = [
            Fq::fq(field_modulus.clone(), BigInt::from(1)),
            Fq::fq(field_modulus.clone(), BigInt::from_str("21888242871839275222246405745257275088696311157297823662689037894645226208582").unwrap())
        ];

        let x1 = Fq::fq(field_modulus.clone(),BigInt::from_str("10857046999023057135944570762232829481370756359578518086990519993285655852781").unwrap());
        let y1 = Fq::fq(field_modulus.clone(),BigInt::from_str("11559732032986387107991004021392285783925812861821192530917403151452391805634").unwrap());
        let x2 = Fq::fq(field_modulus.clone(),BigInt::from_str("8495653923123431417604973247489272438418190587263600148770280649306958101930").unwrap());
        let y2 = Fq::fq(field_modulus.clone(),BigInt::from_str("4082367875863433681332203403145435568316851327593401208105741076214120093531").unwrap());

        let p1 = [
                            Fq2::fq2([x1.clone(), y1.clone()], non_residue.clone(), forbenius_coeffs_c1.clone()),
                            Fq2::fq2([x2.clone(), y2.clone()], non_residue.clone(), forbenius_coeffs_c1.clone()),
                        ];
        let g2 = G2::new(p1[0].clone(), p1[1].clone(), Fq2::one(non_residue.clone(), field_modulus.clone(), forbenius_coeffs_c1.clone()));

        let r1 = get_fq("33");
        let r2 = get_fq("44");

        let gr1 = g2.mul_scalar(r1.clone()).affine();
        let gr2 = g2.mul_scalar(r2.clone()).affine();

        // println!("gr1-{:?} -- {:?} -- {:?}", g2.x.value, g2.y.value, g2.z.value);
        // println!("gr2-{:?} -- {:?} -- {:?}", gr2.x.value, gr2.y.value, gr2.z.value);

        let grsum1 = gr1 + gr2;
        let r1r2 = r1+r2;

        let grsum2 = g2.mul_scalar(r1r2).affine();
        assert_eq!(grsum1, grsum2)
    }
}