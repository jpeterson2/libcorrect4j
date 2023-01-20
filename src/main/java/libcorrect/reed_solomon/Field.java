/*
 * libcorrect4j
 * Field.java
 * Created from include/correct/reed-solomon/field.h @ https://github.com/quiet/libcorrect
 */

package libcorrect.reed_solomon;

import org.checkerframework.checker.signedness.qual.Unsigned;

public class Field {
    private byte[] expTable;
    private byte[] logTable;

    public byte exp(@Unsigned int i) {
        return expTable[i];
    }

    public byte log(@Unsigned int i) {
        return logTable[i];
    }

    public Field(short primitivePoly) {
        // in GF(2^8)
        // log and exp
        // bits are in GF(2), compute alpha^val in GF(2^8)
        // exp should be of size 512 so that it can hold a "wraparound" which prevents some modulo ops
        // log should be of size 256. no wraparound here, the indices into this table are field elements
        expTable = new byte[512];
        logTable = new byte[256];

        // assume alpha is a primitive element, p(x) (primitive_poly) irreducible in GF(2^8)
        // addition is xor
        // subtraction is addition (also xor)
        // e.g. x^5 + x^4 + x^4 + x^2 + 1 = x^5 + x^2 + 1
        // each row of exp contains the field element found by exponentiating
        //   alpha by the row index
        // each row of log contains the coefficients of
        //   alpha^7 + alpha^6 + alpha^5 + alpha^4 + alpha^3 + alpha^2 + alpha + 1
        // as 8 bits packed into one byte

        short element = 1;
        expTable[0] = (byte)element;
        logTable[0] = (byte)0;  // really, it's undefined. we shouldn't ever access this

        for (short i = 1; i < 512; i++) {
            element = (short) (Short.toUnsignedInt(element) * 2);
            element = (short) (Short.toUnsignedInt(element) > 255 ? Short.toUnsignedInt(element) ^ Short.toUnsignedInt(primitivePoly) : Short.toUnsignedInt(element));
            expTable[Short.toUnsignedInt(i)] = (byte) element;
            if (Short.toUnsignedInt(i) < 256) {
                logTable[Short.toUnsignedInt(element)] =(byte) i;
            }
        }
    }

    public byte fieldMulLogElement(byte l, byte r) {
        // like field_mul_log, but returns a field_element_t
        // because we are doing lookup here, we can safely skip the wrapover check
        short res = (short) (Short.toUnsignedInt((short) Byte.toUnsignedInt(l)) + Short.toUnsignedInt((short) Byte.toUnsignedInt(r)));
        return expTable[Short.toUnsignedInt(res)];
    }

    public byte fieldAdd(byte l, byte r) {
        return (byte) (Byte.toUnsignedInt(l) ^ Byte.toUnsignedInt(r));
    }

    public byte fieldSub(byte l, byte r) {
        return (byte) (Byte.toUnsignedInt(l) ^ Byte.toUnsignedInt(r));
    }

    public byte fieldSum(byte elem_U, @Unsigned int n_U) {
        // we'll do a closed-form expression of the sum, although we could also
        //   choose to call field_add n times

        // since the sum is actually the bytewise XOR operator, this suggests two
        // kinds of values: n odd, and n even

        // if you sum once, you have coeff
        // if you sum twice, you have coeff XOR coeff = 0
        // if you sum thrice, you are back at coeff
        // an even number of XORs puts you at 0
        // an odd number of XORs puts you back at your value

        // so, just throw away all the even n
        return (byte) (Integer.remainderUnsigned(n_U, 2) != 0 ? Byte.toUnsignedInt(elem_U) : 0);
    }

    public byte fieldMul(byte l, byte r) {
        if(Byte.toUnsignedInt(l) == 0 || Byte.toUnsignedInt(r) == 0) {
            return 0;
        }
        // multiply two field elements by adding their logarithms.
        // yep, get your slide rules out
        short res = (short)(Short.toUnsignedInt((short)Byte.toUnsignedInt(logTable[Byte.toUnsignedInt(l)])) +
                              Short.toUnsignedInt((short)Byte.toUnsignedInt(logTable[Byte.toUnsignedInt(r)])));

        // if coeff exceeds 255, we would normally have to wrap it back around
        // alpha^255 = 1; alpha^256 = alpha^255 * alpha^1 = alpha^1
        // however, we've constructed exponentiation table so that
        //   we can just directly lookup this result
        // the result must be clamped to [0, 511]
        // the greatest we can see at this step is alpha^255 * alpha^255
        //   = alpha^510
        return expTable[Short.toUnsignedInt(res)];
    }

    public byte fieldDiv(byte l, byte r) {
        if(Byte.toUnsignedInt(l) == 0) {
            return 0;
        }

        if(Byte.toUnsignedInt(r) == 0) {
            // XXX ???
            throw new IllegalArgumentException("fieldDiv: illegal dividor r==0");
        }

        // division as subtraction of logarithms

        // if rcoeff is larger, then log[l] - log[r] wraps under
        // so, instead, always add 255. in some cases, we'll wrap over, but
        // that's ok because the exp table runs up to 511.
        short res = (short)(Short.toUnsignedInt((short)255) +
                              Short.toUnsignedInt((short)Byte.toUnsignedInt(logTable[Byte.toUnsignedInt(l)])) -
                              Short.toUnsignedInt((short)Byte.toUnsignedInt(logTable[Byte.toUnsignedInt(r)])));
        return expTable[Short.toUnsignedInt(res)];
    }

    public byte fieldMulLog(byte l, byte r) {
        // this function performs the equivalent of field_mul on two logarithms
        // we save a little time by skipping the lookup step at the beginning
        short res = (short) (Short.toUnsignedInt((short) Byte.toUnsignedInt(l)) + Short.toUnsignedInt((short) Byte.toUnsignedInt(r)));

        // because we arent using the table, the value we return must be a valid logarithm
        // which we have decided must live in [0, 255] (they are 8-bit values)
        // ensuring this makes it so that multiple muls will not reach past the end of the
        // exp table whenever we finally convert back to an element
        if (Short.toUnsignedInt(res) > 255) {
            return (byte) (Short.toUnsignedInt(res) - 255);
        }
        return (byte) res;
    }

    public byte fieldDivLog(byte l, byte r) {
        // like field_mul_log, this performs field_div without going through a field_element_t
        short res = (short) (Short.toUnsignedInt((short) 255) +
                               Short.toUnsignedInt((short) Byte.toUnsignedInt(l)) -
                               Short.toUnsignedInt((short) Byte.toUnsignedInt(r)));
        if (Short.toUnsignedInt(res) > 255) {
            return (byte) (Short.toUnsignedInt(res) - 255);
        }
        return (byte) res;
    }

    public byte fieldPow(byte elem, @Unsigned int pow) {
        // take the logarithm, multiply, and then "exponentiate"
        // n.b. the exp table only considers powers of alpha, the primitive element
        // but here we have an arbitrary coeff
        byte log_U = this.logTable[Byte.toUnsignedInt(elem)];
        int resLog = Byte.toUnsignedInt(log_U) * pow;
        int mod = resLog % 255; //TODO: unsigned modulus library function I believe
        if (mod < 0) {
            mod += 255;
        }
        return expTable[mod];
    }

}
