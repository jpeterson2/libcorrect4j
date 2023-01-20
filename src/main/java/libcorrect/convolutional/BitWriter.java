/*
 * libcorrect4j
 * BitReader.java
 * Created from src/correct/convolutional/bit.c @ https://github.com/quiet/libcorrect
 */

package libcorrect.convolutional;

import org.checkerframework.checker.signedness.qual.Unsigned;

public class BitWriter {
    private @Unsigned byte currentByte_U;
    private @Unsigned int currentByteLen_U;
    private @Unsigned byte[] bytes_U;
    private @Unsigned int byteIndex_U;
    private @Unsigned long len_U;

    public BitWriter(@Unsigned byte[] bytes_U, @Unsigned long len_U) {
        if (bytes_U != null) {
            reconfigure(bytes_U, len_U);
        }
    }

    public void reconfigure(@Unsigned byte[] bytes_U, @Unsigned long len_U) {
        this.bytes_U = bytes_U;
        this.len_U = len_U;

        currentByte_U = 0;
        currentByteLen_U = 0;
        byteIndex_U = 0;
    }

    public void write(@Unsigned byte val_U, @Unsigned int n_U) {
        for (long j = 0; Long.compareUnsigned(j, Integer.toUnsignedLong(n_U)) < 0; j++) {
            write1(val_U);
            val_U = (byte) (Byte.toUnsignedInt(val_U) >> 1);
        }
    }

    public void write1(@Unsigned byte val_U) {
        currentByte_U = (byte) (Byte.toUnsignedInt((byte) (Byte.toUnsignedInt(currentByte_U) | Byte.toUnsignedInt(val_U) & 1)));
        if (++currentByteLen_U == 8) {
            // 8 bits in a byte -- move to the next byte
            bytes_U[byteIndex_U++] = currentByte_U;
            currentByteLen_U = 0;
            currentByte_U = 0;
        } else {
            currentByte_U = (byte) (Byte.toUnsignedInt(currentByte_U) << 1);
        }
    }

    public void writeBitlist(@Unsigned byte[] l_U, @Unsigned long len_U) {
        int lIndex = 0;
        // first close the current byte
        // we might have been given too few elements to do that. be careful.
        long closeLen_U = Integer.toUnsignedLong(8 - this.currentByteLen_U);
        closeLen_U = Long.compareUnsigned(closeLen_U, len_U) < 0 ? closeLen_U : len_U;

        short b_U = (short) Byte.toUnsignedInt(this.currentByte_U);

        for (long i = 0; Long.compareUnsigned(i, closeLen_U) < 0; i++) {
            b_U = (short) (Short.toUnsignedInt(b_U) | Byte.toUnsignedInt(l_U[lIndex + ((int) i)]));
            b_U = (short) (Short.toUnsignedInt(b_U) << 1);
        }


        lIndex += (int) closeLen_U;
        len_U -= closeLen_U;

        @Unsigned byte[] bytes_U = this.bytes_U;
        int byteIndex_U = this.byteIndex_U;

        if (Integer.toUnsignedLong(this.currentByteLen_U) + closeLen_U == 8) {
            b_U = (short) (Short.toUnsignedInt(b_U) >> 1);
            bytes_U[(int) byteIndex_U] = (byte) b_U;
            byteIndex_U++;
        } else {
            this.currentByte_U = (byte) b_U;
            this.currentByteLen_U = (int) (Integer.toUnsignedLong(this.currentByteLen_U) + closeLen_U);
            return;
        }

        long fullBytes_U = Long.divideUnsigned(len_U, 8);

        for (long i_U = 0; Long.compareUnsigned(i_U, fullBytes_U) < 0; i_U++) {
            bytes_U[(int) byteIndex_U] = (byte) (Byte.toUnsignedInt(l_U[lIndex]) << 7 |
                    Byte.toUnsignedInt(l_U[lIndex + 1]) << 6 |
                    Byte.toUnsignedInt(l_U[lIndex + 2]) << 5 |
                    Byte.toUnsignedInt(l_U[lIndex + 3]) << 4 |
                    Byte.toUnsignedInt(l_U[lIndex + 4]) << 3 |
                    Byte.toUnsignedInt(l_U[lIndex + 5]) << 2 |
                    Byte.toUnsignedInt(l_U[lIndex + 6]) << 1 |
                    Byte.toUnsignedInt(l_U[lIndex + 7]));
            byteIndex_U += 1;
            lIndex += 8;
        }

        len_U -= 8 * fullBytes_U;

        b_U = 0;
        for (long i = 0; Long.compareUnsigned(i, len_U) < 0; i++) {
            b_U = (short) (Short.toUnsignedInt(b_U) | Byte.toUnsignedInt(l_U[lIndex + ((int) i)]));
            b_U = (short) (Short.toUnsignedInt(b_U) << 1);
        }

        this.currentByte_U = (byte) b_U;
        this.byteIndex_U = byteIndex_U;
        this.currentByteLen_U = (int) len_U;
    }

    public void writeBitlistReversed(@Unsigned byte[] l_U, @Unsigned long len_U) {
        int lIndex = 0;
        lIndex += (int) (len_U - 1);

        @Unsigned byte[] bytes_U = this.bytes_U;
        int byteIndex_U = this.byteIndex_U;
        short b_U;

        if (this.currentByteLen_U != 0) {
            long closeLen_U = Integer.toUnsignedLong(8 - this.currentByteLen_U);
            closeLen_U = Long.compareUnsigned(closeLen_U, len_U) < 0 ? closeLen_U : len_U;

            b_U = (short) Byte.toUnsignedInt(this.currentByte_U);

            for (long i = 0; Long.compareUnsigned(i, closeLen_U) < 0; i++) {
                b_U = (short) (Short.toUnsignedInt(b_U) | Byte.toUnsignedInt(l_U[lIndex]));
                b_U = (short) (Short.toUnsignedInt(b_U) << 1);
                lIndex--;
            }

            len_U -= closeLen_U;

            if (Integer.toUnsignedLong(this.currentByteLen_U) + closeLen_U == 8) {
                b_U = (short) (Short.toUnsignedInt(b_U) >> 1);
                bytes_U[(int) byteIndex_U] = (byte) b_U;
                byteIndex_U++;
            } else {
                this.currentByte_U = (byte) b_U;
                this.currentByteLen_U = (int) (Integer.toUnsignedLong(this.currentByteLen_U) + closeLen_U);
                return;
            }
        }

        long fullBytes_U = Long.divideUnsigned(len_U, 8);

        for (long i = 0; Long.compareUnsigned(i, fullBytes_U) < 0; i++) {
            bytes_U[(int) byteIndex_U] = (byte) (Byte.toUnsignedInt(l_U[lIndex]) << 7 |
                    Byte.toUnsignedInt(l_U[lIndex - 1]) << 6 |
                    Byte.toUnsignedInt(l_U[lIndex - 2]) << 5 |
                    Byte.toUnsignedInt(l_U[lIndex - 3]) << 4 |
                    Byte.toUnsignedInt(l_U[lIndex - 4]) << 3 |
                    Byte.toUnsignedInt(l_U[lIndex - 5]) << 2 |
                    Byte.toUnsignedInt(l_U[lIndex - 6]) << 1 |
                    Byte.toUnsignedInt(l_U[lIndex - 7]));
            byteIndex_U += 1;
            lIndex += -8;
        }

        len_U -= 8 * fullBytes_U;

        b_U = 0;
        for (long i = 0; Long.compareUnsigned(i, len_U) < 0; i++) {
            b_U = (short) (Short.toUnsignedInt(b_U) | Byte.toUnsignedInt(l_U[lIndex]));
            b_U = (short) (Short.toUnsignedInt(b_U) << 1);
            lIndex--;
        }

        this.currentByte_U = (byte) b_U;
        this.byteIndex_U = byteIndex_U;
        this.currentByteLen_U = (int) len_U;
    }

    public void flushByte() {
        if (currentByteLen_U != 0) {
            currentByte_U = (byte) (Byte.toUnsignedInt(this.currentByte_U) << 8 - this.currentByteLen_U);
            bytes_U[byteIndex_U] = currentByte_U;
            byteIndex_U++;
            currentByteLen_U = 0;
        }
    }

    public @Unsigned int length() {
        return byteIndex_U;
    }
}
