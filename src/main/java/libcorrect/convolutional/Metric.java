/*
 * libcorrect4j
 * Convolutional.java
 * Created from     include/correct/convolutional/metric.h
 *                  src/correct/convolutional/metric.c @ https://github.com/quiet/libcorrect
 */
package libcorrect.convolutional;

import static java.lang.Integer.bitCount;
import org.checkerframework.checker.signedness.qual.Unsigned;

public class Metric {
    /**
     * implemented as population count of x XOR y
     */
    public static @Unsigned short distance(@Unsigned int x_U, @Unsigned int y_U) {
        return (short)bitCount(x_U ^ y_U);
    }

    public static @Unsigned short softDistanceLinear(@Unsigned int hardX_U, @Unsigned byte[] softY_U, @Unsigned long len_U, @Unsigned int shift) {
        short dist_U = 0;
        for(int i = 0; Long.compareUnsigned(Integer.toUnsignedLong(i), len_U) < 0; i++) {
            int softX_U = (byte)0 - (hardX_U & 1) & 0xff;
            hardX_U >>>= 1;
            int d = Byte.toUnsignedInt(softY_U[i+shift]) - softX_U;
            dist_U = (short)(Short.toUnsignedInt(dist_U) + (d < 0 ? -d : d));
        }
        return dist_U;
    }
    /**
     * since euclidean dist is sqrt(a^2 + b^2 + ... + n^2), the square is just
     * a^2 + b^2 + ... + n^2
     */
    public static @Unsigned short softDistanceQuadratic(@Unsigned int hardX_U, @Unsigned byte[] softY_U, @Unsigned long len_U, @Unsigned int shift) {
        short dist_U = 0;
        for(int i = 0; Long.compareUnsigned(Integer.toUnsignedLong(i), len_U) < 0; i++) {
            // first, convert hard_x to a soft measurement (0 -> 0, 1 - > 255)
            int softX_U = (hardX_U & 1) != 0 ? 255 : 0;
            hardX_U >>>= 1;
            int d = Byte.toUnsignedInt(softY_U[i+shift]) - softX_U;
            dist_U = (short)(Short.toUnsignedInt(dist_U) + d * d);
        }
        return (short)(Short.toUnsignedInt(dist_U) >>> 3);
    }

}
