use std::ops;
use num_bigint::BigInt;
use crate::field::fq::Fq;
use crate::field::fq2::Fq2;


#[derive(Debug, Clone)]
pub struct G1Affine {
    x: Fq,
    y: Fq
}

// has three points x,y,z since we use jacobian coordinates
// X = x/z^2, Y = y/z^3
#[derive(Debug, Clone)]
pub struct G1 {
    pub x: Fq,
    pub y: Fq,
    pub z: Fq
}

impl G1 {
    pub fn new(x: Fq, y: Fq, z:Fq) -> Self {
        Self {
            x, y, z
        }
    }

    pub fn double(&self) -> Self {
        if self.z.value == BigInt::from(0) {
            return self.clone();
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
        return G1::new(x3, y3, z3);
    }

    pub fn mul_scalar(&self, base: Fq) -> Self {
        let d = base.clone();
        let mut q = G1::new(
            Fq::zero(d.field_modulus.clone()),
            Fq::zero(d.field_modulus.clone()),
            Fq::zero(d.field_modulus.clone())
        );
        let r = self.clone();
        let mut found_one = false;

        let d_bit_len = d.value.bits();
        for i in (0..d_bit_len).rev() {
            if found_one{
                q = q.double();
            }
            if d.value.bit(i) {
                found_one = true;
                q = q + r.clone();
            }
        }
        return q;
    }

    pub fn affine(&self, fq_field_modulus: BigInt) -> G1Affine {
        if self.z.value == BigInt::from(0) {
            return G1Affine{x: self.x.clone(), y: self.y.clone()}
        }
        let z_inv = Fq::one(fq_field_modulus)/self.z.clone();
        let z_inv2 = z_inv.clone() * z_inv.clone();
        let x = self.x.clone() * z_inv2.clone();

        let z_inv3 = z_inv2.clone() * z_inv.clone();
        let y = self.y.clone() * z_inv3.clone();

        return G1Affine{x, y};
    }
}

impl ops::Add<G1> for G1 {
    type Output = G1;

    fn add(self, rhs: G1) -> Self::Output {
        if self.z.value == BigInt::from(0) {
            return rhs.clone()
        }
        if rhs.z.value == BigInt::from(0) {
            return self.clone()
        }
        let x1 = self.x.clone();
        let y1 = self.y.clone();
        let z1 = self.z.clone();
        let x2 = rhs.x.clone();
        let y2 = rhs.y.clone();
        let z2 = rhs.z.clone();

        let z1_2 = z1.clone() * z1.clone();
        let z2_2 = z2.clone() * z2.clone();

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

        return G1::new(x3, y3, z3);
    }
}

impl PartialEq for G1 {
    fn eq(&self, other: &Self) -> bool {
        if self.z.value == BigInt::from(0) {
            return other.z.value == BigInt::from(0)
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

impl ops::Neg for G1 {
    type Output = G1;

    fn neg(self) -> Self::Output {
        return G1::new(self.x, -self.y, self.z)
    }
}

impl ops::Sub<G1> for G1 {
    type Output = G1;

    fn sub(self, rhs: G1) -> Self::Output {
        return self.clone() + (-rhs.clone())
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
    fn test_g1_add() {
        let fq_field_modulus = BigInt::from_str("21888242871839275222246405745257275088696311157297823662689037894645226208583").unwrap();
        let mut gr1 = G1::new(Fq::one(fq_field_modulus.clone()), get_fq("2"), Fq::one(fq_field_modulus.clone()));
        let mut gr2 = G1::new(Fq::one(fq_field_modulus.clone()), get_fq("2"), Fq::one(fq_field_modulus.clone()));
        let x = G1::new(Fq::one(fq_field_modulus.clone()), get_fq("2"), Fq::one(fq_field_modulus.clone()));

        let r1 = get_fq("33");
        let r2 = get_fq("44");

        let gr1 = gr1.mul_scalar(r1.clone());

        let gr2 = gr2.mul_scalar(r2.clone());

        let r1r2 = r1.clone() + r2.clone();

        let grsum1 = gr1.clone() + gr2.clone();
        let grsum2 = x.mul_scalar(r1r2.clone());

        assert_eq!(grsum1, grsum2);

        let a = grsum1.affine(fq_field_modulus.clone());

        assert_eq!(a.x.value, BigInt::from_str("21526464323725832663882905544083280657770325585151797133383551854196089356032").unwrap());
        assert_eq!(a.y.value, BigInt::from_str("8545759555567142326482563981456384114560528812235279370284511019768507753138").unwrap());
    }
}
