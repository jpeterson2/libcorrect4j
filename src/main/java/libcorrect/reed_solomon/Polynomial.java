/*
 * libcorrect4j
 * Polynomial.java
 * Created from src/correct/reed-solomon/polynomial.c @ https://github.com/quiet/libcorrect
 */

package libcorrect.reed_solomon;

import java.util.Arrays;
import org.checkerframework.checker.signedness.qual.Unsigned;

public class Polynomial implements Cloneable {
    private byte[] coeff;
    private int order;

    /**
     * Create uninitialized polyno,ial
     * @param order
     */
    public Polynomial(int order) {
        coeff = new byte[order + 1]; // makes this compile-time rather than runtime (NegativeArraySizeException)
        this.order = order;
    }

    /**
     * Create polynomial from coefficients
     * @param order
     * @param coeff
     */
    public Polynomial(int order, byte[] coeff) {
        this.order = order;
        this.coeff = Arrays.copyOf(coeff, this.order +1);
    }

    /**
     * Create polynomial from roots
     * @param field
     * @param nRoots
     * @param roots
     */
    public Polynomial(Field field, int nRoots, byte[] roots) {
        Polynomial l = new Polynomial(1);

        Polynomial[] r = new Polynomial[2];
        // we'll keep two temporary stores of rightside polynomial
        // each time through the loop, we take the previous result and use it as new rightside
        // swap back and forth (prevents the need for a copy)

        r[0] = new Polynomial(nRoots);
        r[1] = new Polynomial(nRoots);
        r[0].order = 1;

        int rCoeffRes = 0;

        // initialize the result with x + roots[0]
        r[rCoeffRes].coeff[0] = roots[0];
        r[rCoeffRes].coeff[1] = (byte) 1;

        // initialize lcoeff[1] with x
        // we'll fill in the 0th order term in each loop iter
        l.coeff[1] = (byte) 1;

        // loop through, using previous run's result as the new right hand side
        // this allows us to multiply one group at a time
        for (int i = 1; Integer.compare(i, nRoots) < 0; i++) {
            l.coeff[0] =  roots[i];
            int nextRCoeff = rCoeffRes;
            rCoeffRes = (rCoeffRes + 1) % 2;
            r[rCoeffRes].order = i + 1;
            mul(field, l, r[nextRCoeff], r[rCoeffRes]);
        }

        coeff = Arrays.copyOf(r[rCoeffRes].coeff, nRoots + 1);
        order = nRoots;

    }


    public void copyCoeff(Polynomial p) {
        System.arraycopy(p.coeff, 0, coeff, 0, p.coeff.length);
    }

    public int getOrder() {
        return order;
    }

    public void setOrder(int v) {
        order = v;
    }

    public byte[] getCoeff() {
        return coeff;
    }

    public byte getCoeff(int i) {
        return coeff[i];
    }

    public void setCoeff(int i, byte v) {
        coeff[i] = v;
    }

    public void flushCoeff() {
        Arrays.fill(coeff, (byte)0); // I don't think actually a bug but there's maybe a way to annotate this so it is understood
    }



    /**
     * Perform an element-wise multiplication of two polynomials
     * if you want a full multiplication, then make res.order = l.order + r.order
     * but if you just care about a lower order, e.g. mul mod x^i, then you can select
     * fewer coefficients
     * @param field
     * @param l
     * @param r
     * @param res
     */

    public static void mul(Field field, Polynomial l, Polynomial r, Polynomial res) {
         Arrays.fill(res.coeff, (byte) 0);
        for (int i = 0; Integer.compare(i, l.order) <= 0; i++) {
            if (Integer.compare(i, res.order) > 0) {
                continue;
            }
            int jLimit_U = Integer.compare(r.order, res.order - i) > 0 ? res.order - i : r.order;
            for (int j = 0; Integer.compare(j, jLimit_U) <= 0; j++) {
                // e.g. alpha^5*x * alpha^37*x^2 --> alpha^42*x^3
                res.coeff[i + j] = field.fieldAdd(res.coeff[i + j],
                        field.fieldMul(l.coeff[i], r.coeff[j]));
            }
        }
    }

    /**
     * Find the polynomial remainder of dividend mod divisor
     * do long division and return just the remainder (written to mod)
     * @param field
     * @param dividend
     * @param divisor
     * @param mod
     */
    public static void mod(Field field, Polynomial dividend, Polynomial divisor, Polynomial mod) {

        if (Integer.compare(mod.order, dividend.order) < 0) {
            // mod.order must be >= dividend.order (scratch space needed)
            // this is an error -- catch it in debug?
            return;
        }
        // initialize remainder as dividend
        mod.coeff = Arrays.copyOf(dividend.coeff, dividend.order + 1);


        // XXX make sure divisor[divisor_order] is nonzero
        byte divisorLeading = field.log(Byte.toUnsignedInt(divisor.coeff[divisor.order]));

        // long division steps along one order at a time, starting at the highest order
        for (int i = dividend.order; Integer.compare(i, 0) > 0; i--) {
            // look at the leading coefficient of dividend and divisor
            // if leading coefficient of dividend / leading coefficient of divisor is q
            //   then the next row of subtraction will be q * divisor
            // if order of q < 0 then what we have is the remainder and we are done
            if (Integer.compare(i, divisor.order) < 0) {
                break;
            }
            if (Byte.toUnsignedInt(mod.coeff[i]) == 0) {
                continue;
            }
            int qOrder = i - divisor.order;
            byte qCoeff = field.fieldDivLog(field.log(Byte.toUnsignedInt(mod.coeff[i])), divisorLeading);

            // now that we've chosen q, multiply the divisor by q and subtract from
            //   our remainder. subtracting in GF(2^8) is XOR, just like addition
            for (int j = 0; Integer.compare(j, divisor.order) <= 0; j++) {
                if (Byte.toUnsignedInt(divisor.coeff[j]) == 0) {
                    continue;
                }
                // all of the multiplication is shifted up by q_order places
                mod.coeff[j + qOrder] = field.fieldAdd(mod.coeff[j + qOrder],
                        field.fieldMulLogElement(field.log(Byte.toUnsignedInt(divisor.coeff[j])), qCoeff));

            }
        }
    }

    /**
     * if f(x) = a(n)*x^n + ... + a(1)*x + a(0)
     * then f'(x) = n*a(n)*x^(n-1) + ... + 2*a(2)*x + a(1)
     * where n*a(n) = sum(k=1, n, a(n)) e.g. the nth sum of a(n) in GF(2^8)
     * @param field
     * @param poly
     * @param der
     */
    public static void formalDerivative(Field field, Polynomial poly, Polynomial der) {

        // assumes der.order = poly.order - 1
        Arrays.fill(der.coeff, (byte) 0);
        for (int i = 0; Integer.compare(i, der.order) <= 0; i++) {
            // we're filling in the ith power of der, so we look ahead one power in poly
            // f(x) = a(i + 1)*x^(i + 1) -> f'(x) = (i + 1)*a(i + 1)*x^i
            // where (i + 1)*a(i + 1) is the sum of a(i + 1) (i + 1) times, not the product
            der.coeff[i] = field.fieldSum(poly.coeff[i + 1], i + 1);
        }
    }

    /**
     * Evaluate the polynomial poly at a particular element val
     * @param field
     * @param poly
     * @param val
     * @return
     */
    public static byte eval(Field field, Polynomial poly, byte val) {

        if (Byte.toUnsignedInt(val) == 0) {
            return poly.coeff[0];
        }

        byte res = 0;

        // we're going to start at 0th order and multiply by val each time
        byte valExponentiated = field.log(1);
        byte valLog = field.log(Byte.toUnsignedInt(val));

        for (int i = 0; Integer.compare(i, poly.order) <= 0; i++) {
            if (Byte.toUnsignedInt(poly.coeff[i]) != 0) {
                // multiply-accumulate by the next coeff times the next power of val
                res = field.fieldAdd(res, field.fieldMulLogElement(field.log(Byte.toUnsignedInt(poly.coeff[i])),
                                                                              valExponentiated));
            }
            // now advance to the next power
            valExponentiated = field.fieldMulLog(valExponentiated, valLog);
        }
        return res;
    }

    /**
     * Evaluate the polynomial poly at a particular element val
     * in this case, all of the logarithms of the successive powers of val have been precalculated
     * this removes the extra work we'd have to do to calculate val_exponentiated each time
     * if this function is to be called on the same val multiple times
     * @param field
     * @param valExp
     * @return
     */
    public byte evalLut(Field field, byte[] valExp) {
        if (Byte.toUnsignedInt(valExp[0]) == 0) {
            return coeff[0];
        }

        byte res = 0;

        for (int i = 0; Integer.compare(i, order) <= 0; i++) {
            if (Byte.toUnsignedInt(coeff[i]) != 0) {
                // multiply-accumulate by the next coeff times the next power of val
                res = field.fieldAdd(res, field.fieldMulLogElement(field.log(Byte.toUnsignedInt(coeff[i])), valExp[i]));
            }
        }
        return res;
    }

    /**
     * Evaluate the log_polynomial poly at a particular element val
     * like polynomial_eval_lut, the logarithms of the successive powers of val have been
     * precomputed
     * @param field
     * @param valExp
     * @return
     */
    public byte evalLogLut(Field field, byte[] valExp) {
        if (Byte.toUnsignedInt(valExp[0]) == 0) {
            if (Byte.toUnsignedInt(coeff[0]) == 0) {
                // special case for the non-existant log case
                return 0;
            }
            return field.exp(Byte.toUnsignedInt(coeff[0]));
        }

        byte res = 0;

        for (int i = 0; Integer.compare(i, order) <= 0; i++) {
            // using 0 as a sentinel value in log -- log(0) is really -inf
            if (Byte.toUnsignedInt(coeff[i]) != 0) {
                // multiply-accumulate by the next coeff times the next power of val
                res = field.fieldAdd(res, field.fieldMulLogElement(coeff[i], valExp[i]));
            }
        }
        return res;
    }

    /**
     * Create the lookup table of successive powers of val used by polynomial_eval_lut
     * @param field
     * @param val
     * @param order
     * @param valExp
     */
    public static void buildExpLut(Field field, byte val, int order, byte[] valExp) {
        byte valExponentiated = field.log(1);
        byte valLog = field.log(Byte.toUnsignedInt(val));
        for (int i = 0; Integer.compare(i, order) <= 0; i++) {
            if (Byte.toUnsignedInt(val) == 0) {
                valExp[i] = (byte) 0;
            } else {
                valExp[i] = valExponentiated;
                valExponentiated = field.fieldMulLog(valExponentiated, valLog);
            }
        }
    }

    /**
     * Initialize polynomial from roots
     * @param field
     * @param nRoots
     * @param roots
     * @param scratch
     */
    public void initFromRoots(Field field, int nRoots,  byte[] roots,  Polynomial[] scratch) {
        Polynomial l = new Polynomial(1);

        // we'll keep two temporary stores of rightside polynomial
        // each time through the loop, we take the previous result and use it as new rightside
        // swap back and forth (prevents the need for a copy)
        Polynomial[] r = new Polynomial[2];
        r[0] = scratch[0];  // clone ?
        r[1] = scratch[1];  // clone ?
        int rCoeffRes = 0;

        // initialize the result with x + roots[0]
        r[rCoeffRes].coeff[1] = (byte) 1;
        r[rCoeffRes].coeff[0] = roots[0];
        r[rCoeffRes].order = 1;

        // initialize lcoeff[1] with x
        // we'll fill in the 0th order term in each loop iter
        l.coeff[1] = (byte) 1;

        // loop through, using previous run's result as the new right hand side
        // this allows us to multiply one group at a time
        for (int i = 1; Integer.compare(i, nRoots) < 0; i++) {
            l.coeff[0] = roots[i];
            int nextRCoeff = rCoeffRes;
            rCoeffRes = (rCoeffRes + 1) % 2;
            r[rCoeffRes].order = i + 1;
            mul(field, l, r[nextRCoeff], r[rCoeffRes]);
        }

        coeff = Arrays.copyOf(r[rCoeffRes].coeff, nRoots + 1);
        order = nRoots;
    }

}

