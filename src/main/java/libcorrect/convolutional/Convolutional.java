/*
 * libcorrect4j
 * Convolutional.java
 * Created from src/correct/convolutional/convolutional.c
 * 				include/correct/convolutional.h
 * 				src/correct/convolutional/encode.c
 * 				src/correct/convolutional/decode.c  @ https://github.com/quiet/libcorrect
 */

package libcorrect.convolutional;


import java.util.Arrays;

import org.checkerframework.checker.signedness.qual.Unsigned;

import static libcorrect.convolutional.Metric.distance;

public class Convolutional {
    // Convolutional Codes
    // Convolutional polynomials are 16 bits wide
    public final static short[] correctConvR126Polynomial = {073, 061};
    public final static short[] correctConvR127Polynomial = {0161, 0127};
    public final static short[] correctConvR128Polynomial = {0225, 0373};
    public final static short[] correctConvR129Polynomial = {0767, 0545};
    public final static short[] correctConvR136Polynomial = {053, 075, 047};
    public final static short[] correctConvR137Polynomial = {0137, 0153, 0121};
    public final static short[] correctConvR138Polynomial = {0333, 0257, 0351};
    public final static short[] correctConvR139Polynomial = {0417, 0627, 0675};

    public final static int CORRECT_SOFT_LINEAR = 0;
    public final static int CORRECT_SOFT_QUADRATIC = CORRECT_SOFT_LINEAR + 1;

    // Maximum of unsigned integral types.
    private final static int UINT8_MAX = 255;
    private final static int UINT16_MAX = 65_535;
    private final static byte SOFT_MAX = (byte) UINT8_MAX;
    private final static short DISTANCE_MAX = (short) UINT16_MAX;


    private final @Unsigned int[] table_U;                /* size 2**order */
    private final @Unsigned int rate_U;                   /* e.g. 2, 3...  */
    private final @Unsigned int order_U;                  /* e.g. 7, 9...	 */
    private final @Unsigned int numstates_U;              /* 2**order 	 */
    private final BitWriter bitWriter;
    private final BitReader bitReader;

    private boolean hasInitDecode;
    private @Unsigned short[] distances_U;
    private PairLookup pairLookup;
    private int softMeasurement;
    private HistoryBuffer historyBuffer;
    private ErrorBuffer errorBuffer;

    /**
     * Convolutional encoder/decoder constructor
     * Allocates and initializes an encoder/decoder for
     * a convolutional code with the given parameters. This function expects that
     * poly will contain inv_rate elements. E.g., to create a conv. code instance
     * with rate 1/2, order 7 and polynomials 0161, 0127, call
     * new Convolutional (2, 7, {0161, 0127});
     *
     * @param r    inverted rate
     * @param o    order
     * @param p    polynomials
     * @throws IllegalArgumentException if requested encoder/decoder cannot be created
     */
    public Convolutional(@Unsigned int r, @Unsigned int o, @Unsigned short[] p) throws IllegalArgumentException {
        if (Integer.compareUnsigned(o, Integer.SIZE) > 0) {
            // XXX turn this into an error code
            // printf("order must be smaller than 8 * sizeof(shift_register_t)\n");
            throw new IllegalArgumentException("Convolutional: order must be smaller than 8 * sizeof(shift_register_t)");
        }
        if (Long.compareUnsigned(r, 2) < 0) {
            // XXX turn this into an error code
            // printf("rate must be 2 or greater\n");
            throw new IllegalArgumentException("Convolutional: rate must be 2 or greater");
        }
        order_U = o;
        rate_U = r;
        numstates_U = (1 << o);

        table_U = new @Unsigned int[1 << o];
        fillTable(rate_U, order_U, p);

        bitWriter = new BitWriter(null, 0);
        bitReader = new BitReader(null, 0);

        hasInitDecode = false;
    }

    /**
     * encodeLen returns the number of *bits*
     * in a msgLen of given size, in *bytes*. In order to convert
     * this returned length to bytes, save the result of the length
     * modulo 8. If it's nonzero, then the length in bytes is
     * length/8 + 1. If it is zero, then the length is just
     * length/8.
     * @param msgLen        the message length
     * @return              the number of *bits* in *bytes*
     */
    public @Unsigned long encodeLen(@Unsigned long msgLen) {
        long msgbits = 8 * msgLen;
        long encodedbits = rate_U * (msgbits + order_U + 1);
        return encodedbits;
    }

    /**
     * Encode a block of data
     * @param msg   a message to encode
     * @return      encoded message
     */

    // assume that encoded length is long enough?
    public @Unsigned byte[] encode(@Unsigned byte[] msg) {
        int msgLen = msg.length;

        // convolutional code convolves filter coefficients, given by
        //     the polynomial, with some history from our message.
        //     the history is stored as single subsequent bits in shiftregister

        // shift in most significant bit every time, one byte at a time
        // shift register takes most recent bit on right, shifts left
        // poly is written in same order, just & mask message w/ poly
        int shiftregister_U = 0;

        // shiftmask is the shiftregister bit mask that removes bits
        //      that extend beyond order
        // e.g. if order is 7, then remove the 8th bit and beyond
        int shiftmask_U = (1 << order_U) - 1;

        long encodedLenBits = encodeLen(msgLen);
        int encodedLen = (int) (Long.remainderUnsigned(encodedLenBits, 8) != 0 ?
                        Long.divideUnsigned(encodedLenBits, 8) + 1 :
                        Long.divideUnsigned(encodedLenBits, 8));
        @Unsigned byte[] encoded = new @Unsigned byte[encodedLen];
        bitWriter.reconfigure(encoded, encodedLen);
        bitReader.reconfigure(msg, msgLen);

        for (long i = 0; Long.compareUnsigned(i, 8 * msgLen) < 0; i++) {
            // shiftregister has oldest bits on left, newest on right
            shiftregister_U <<= 1;
            shiftregister_U |= bitReader.read(1);
            shiftregister_U &= shiftmask_U;
            // shift most significant bit from byte and move down one bit at a time

            // we do direct lookup of our convolutional output here
            // all of the bits from this convolution are stored in this row
            int out_U = table_U[shiftregister_U];
            bitWriter.write((byte) out_U, rate_U);
        }
        // now flush the shiftregister
        // this is simply running the loop as above but without any new inputs
        // or rather, the new input string is all 0s
        for (long i = 0; Long.compareUnsigned(i, order_U + 1) < 0; i++) {
            shiftregister_U <<= 1;
            shiftregister_U &= shiftmask_U;
            int out_U = table_U[shiftregister_U];
            bitWriter.write((byte) out_U, rate_U);
        }

        // 0-fill any remaining bits on our final byte
        bitWriter.flushByte();

        return encoded;
    }

    /**
     * decodeSoft decodes a block encoded by encode and
     * then modulated/demodulated to 8-bit symbols. This function expects
     * that 1 is mapped to 255 and 0 to 0. An erased symbol should be
     * set to 128. The decoded message may contain errors.
     * @param encoded
     * @param numEncodedBits  should contain the length of encoded in *bits*.
     *       This value need not be an exact multiple of 8. However,
     *       it must be a multiple of the inv_rate used to create
     *       the conv instance.
     * @return  decoded message
     * @throws IllegalArgumentException
     */

    public @Unsigned byte[] decodeSoft(@Unsigned byte[] encoded, @Unsigned long numEncodedBits) throws IllegalArgumentException {
        if (Long.remainderUnsigned(numEncodedBits, this.rate_U) != 0) {
            // XXX turn this into an error code
            // printf("encoded length of message must be a multiple of rate\n");
            throw new IllegalArgumentException("correctConvolutionalDecode: encoded length of message must be a multiple of rate");
        }

        long numEncodedBytes = Long.remainderUnsigned(numEncodedBits, 8) != 0 ? Long.divideUnsigned(numEncodedBits, 8) + 1 :
                Long.divideUnsigned(numEncodedBits, 8);

        return _decode(numEncodedBits, numEncodedBytes,encoded);
    }

     /**
     * @param numEncodedBits  should contain the length of encoded in *bits*.
     *       This value need not be an exact multiple of 8. However,
     *       it must be a multiple of the inv_rate used to create
     *       the conv instance.
     * @return  decoded message
     * @throws IllegalArgumentException
     */
    public @Unsigned byte[] decode(@Unsigned byte[] encoded, @Unsigned long numEncodedBits) throws IllegalArgumentException {
        if (Long.remainderUnsigned(numEncodedBits, this.rate_U) != 0) {
            // XXX turn this into an error code
            // printf("encoded length of message must be a multiple of rate\n");
            throw new IllegalArgumentException("correctConvolutionalDecode: encoded length of message must be a multiple of rate");
        }

        long numEncodedBytes = Long.remainderUnsigned(numEncodedBits, 8) != 0 ? Long.divideUnsigned(numEncodedBits, 8) + 1 :
                Long.divideUnsigned(numEncodedBits, 8);
        bitReader.reconfigure(encoded, numEncodedBytes);

        return _decode(numEncodedBits, numEncodedBytes, null);
    }


    public void decodePrintIter(@Unsigned int iter, @Unsigned int winnerIndex) {
        if (Integer.compareUnsigned(iter, 2_220) < 0) {
            return;
        }
        System.out.println("iteration: " + iter);
        System.out.println("errors:");
        for (int i = 0; Integer.compareUnsigned(i, Integer.divideUnsigned(numstates_U, 2)) < 0; i++) {
            System.out.printf("%2d: %d\n", i, Short.toUnsignedInt(errorBuffer.getWriteError(i)));
        }
        System.out.println();
        System.out.println("history:");
        for (int i = 0; Integer.compareUnsigned(i, Integer.divideUnsigned(numstates_U, 2)) < 0; i++) {
            System.out.printf("%2d: ", i);
            for (int j = 0; Integer.compareUnsigned(j, winnerIndex) <= 0; j++) {
                System.out.printf("%d", historyBuffer.getHistory(j, i) != 0 ? 1 : 0);
            }
            System.out.println();
        }
        System.out.println();
    }

    private @Unsigned byte[] _decode(@Unsigned long numEncodedBits_U, @Unsigned long numEncodedBytes_U, @Unsigned byte[] softEncoded_U) {
        if (!hasInitDecode) {
            long maxErrorPerInput_U = rate_U * Byte.toUnsignedLong(SOFT_MAX);
            int renormalizeInterval_U = (int) Long.divideUnsigned(Short.toUnsignedLong(DISTANCE_MAX), maxErrorPerInput_U);
            decodeInit(5 * order_U, 15 * order_U, renormalizeInterval_U);
        }

        int sets_U = (int) Long.divideUnsigned(numEncodedBits_U, rate_U);
        // XXX fix this vvvvvv
        @Unsigned byte[] msg = new @Unsigned byte[(int) numEncodedBytes_U];
        bitWriter.reconfigure(msg, numEncodedBytes_U);
        errorBuffer.reset();
        historyBuffer.reset();

        // no outputs are generated during warmup
        decodeWarmup(sets_U, softEncoded_U);
        decodeInner(sets_U, softEncoded_U);
        decodeTail(sets_U, softEncoded_U);

        historyBuffer.flush(bitWriter);

        return Arrays.copyOf(msg, bitWriter.length());
    }

    private void decodeInit(@Unsigned int minTraceback_U, @Unsigned int tracebackLength_U, @Unsigned int renormalizeInterval_U) {
        hasInitDecode = true;
        distances_U = new short[1 << rate_U];
        pairLookup = new PairLookup(rate_U, order_U, table_U);
        softMeasurement = CORRECT_SOFT_LINEAR;

        // we limit history to go back as far as 5 * the order of our polynomial
        historyBuffer = new HistoryBuffer(minTraceback_U, tracebackLength_U, renormalizeInterval_U,
                numstates_U / 2, 1 << (order_U - 1)); // TODO: unsigned division?
        errorBuffer = new ErrorBuffer(numstates_U);
    }

    private void decodeWarmup(@Unsigned int sets_U, @Unsigned byte[] soft_U) {
        // first phase: load shiftregister up from 0 (order goes from 1 to conv->order)
        // we are building up error metrics for the first order bits
        for (int i = 0; Long.compareUnsigned(Integer.toUnsignedLong(i), this.order_U - 1) < 0 && Integer.compareUnsigned(i, sets_U) < 0; i++) {
            // peel off rate bits from encoded to recover the same `out` as in the encoding process
            // the difference being that this `out` will have the channel noise/errors applied
            int out_U = 0;
            if (soft_U == null) {
                out_U = bitReader.read(this.rate_U);
            }
            // walk all of the state we have so far
            for (int j = 0; Integer.compareUnsigned(j, 1 << i + 1) < 0; j += 1) {
                int last_U = (int) (j >>> 1);
                short dist_U;

                if (soft_U != null) {
                    if (softMeasurement == CORRECT_SOFT_LINEAR) {
                        dist_U = Metric.softDistanceLinear(table_U[j], soft_U, this.rate_U, i * this.rate_U);
                    } else {
                        dist_U = Metric.softDistanceQuadratic(table_U[j], soft_U, this.rate_U, i * this.rate_U);
                    }
                } else {
                    dist_U = distance(table_U[j], out_U);
                }
                errorBuffer.setWriteError(j, (short) (dist_U + errorBuffer.getReadError(last_U)));
            }
            errorBuffer.swap();
        }
    }

    private void decodeInner(@Unsigned int sets_U, @Unsigned byte[] soft_U) {
        int highbit_U = 1 << order_U - 1;
        for (int i = order_U - 1; Long.compareUnsigned(Integer.toUnsignedLong(i), Integer.toUnsignedLong(sets_U) - order_U + 1) < 0; i++) {
            // lasterrors are the aggregate bit errors for the states of shiftregister for the previous
            // time slice
            if (soft_U != null) {
                if (softMeasurement == CORRECT_SOFT_LINEAR) {
                    for (int j = 0; Integer.compareUnsigned(j, 1 << rate_U) < 0; j++) {
                        distances_U[j] = Metric.softDistanceLinear(j, soft_U, rate_U, i * rate_U);
                    }
                } else {
                    for (int j = 0; Integer.compareUnsigned(j, 1 << rate_U) < 0; j++) {
                        distances_U[j] = (short) Metric.softDistanceQuadratic(j, soft_U, rate_U, i * rate_U);
                    }
                }
            } else {
                int out_U = bitReader.read(rate_U);
                for (int i2 = 0; Integer.compareUnsigned(i2, 1 << rate_U) < 0; i2++) {
                    distances_U[i2] = distance(i2, out_U);
                }
            }

            pairLookup.fillDistance(distances_U);

            // a mask to get the high order bit from the shift register
            int numIter_U = highbit_U << 1;
            // aggregate bit errors for this time slice

            @Unsigned byte[] history_U = historyBuffer.getSlice();
            // walk through all states, ignoring oldest bit
            // we will track a best register state (path) and the number of bit errors at that path at
            // this time slice
            // this loop considers two paths per iteration (high order bit set, clear)
            // so, it only runs numstates/2 iterations
            // we'll update the history for every state and find the path with the least aggregated bit
            // errors

            // now run the main loop
            // we calculate 2 sets of 2 register states here (4 states per iter)
            // this creates 2 sets which share a predecessor, and 2 sets which share a successor
            //
            // the first set definition is the two states that are the same except for the least order
            // bit
            // these two share a predecessor because their high n - 1 bits are the same (differ only by
            // newest bit)
            //
            // the second set definition is the two states that are the same except for the high order
            // bit
            // these two share a successor because the oldest high order bit will be shifted out, and
            // the other bits will be present in the successor
            //
            int highbase_U = highbit_U >>> 1;
            for (int low_U = 0, high_U = highbit_U, base_U = 0; Integer.compareUnsigned(high_U, numIter_U) < 0; low_U += 8, high_U += 8, base_U += 4) {
                // shifted-right ancestors
                // low and low_plus_one share low_past_error
                //   note that they are the same when shifted right by 1
                // same goes for high and high_plus_one
                for (int offset_U = 0, baseOffset_U = 0; Integer.compareUnsigned(baseOffset_U, 4) < 0; offset_U += 2, baseOffset_U += 1) {
                    int lowKey_U = pairLookup.getKey(base_U + baseOffset_U);
                    int highKey_U = pairLookup.getKey(highbase_U + base_U + baseOffset_U);
                    int lowConcatDist_U = pairLookup.getDistance(lowKey_U);
                    int highConcatDist_U = pairLookup.getDistance(highKey_U);

                    short lowPastError_U = errorBuffer.getReadError(base_U + baseOffset_U);
                    short highPastError_U = errorBuffer.getReadError(highbase_U + base_U + baseOffset_U);

                    short lowError_U = (short) ((lowConcatDist_U & 0xffff) + Short.toUnsignedInt(lowPastError_U));
                    short highError_U = (short) ((highConcatDist_U & 0xffff) + Short.toUnsignedInt(highPastError_U));

                    int successor_U = low_U + offset_U;
                    short error_U;
                    byte historyMask_U;
                    if (Short.toUnsignedInt(lowError_U) <= Short.toUnsignedInt(highError_U)) {
                        error_U = lowError_U;
                        historyMask_U = 0;
                    } else {
                        error_U = highError_U;
                        historyMask_U = 1;
                    }
                    errorBuffer.setWriteError(successor_U, error_U);
                    history_U[successor_U] = historyMask_U;

                    int lowPlusOne_U = low_U + offset_U + 1;
                    short lowPlusOneError_U = (short) ((lowConcatDist_U >>> 16) + Short.toUnsignedInt(lowPastError_U));
                    short highPlusOneError_U = (short) ((highConcatDist_U >>> 16) + Short.toUnsignedInt(highPastError_U));

                    short plusOneError_U;
                    int plusOneSuccessor_U = lowPlusOne_U;
                    byte plusOneHistoryMask_U;
                    if (Short.toUnsignedInt(lowPlusOneError_U) <= Short.toUnsignedInt(highPlusOneError_U)) {
                        plusOneError_U = lowPlusOneError_U;
                        plusOneHistoryMask_U = 0;
                    } else {
                        plusOneError_U = highPlusOneError_U;
                        plusOneHistoryMask_U = 1;
                    }

                    errorBuffer.setWriteError(plusOneSuccessor_U, plusOneError_U);
                    history_U[plusOneSuccessor_U] = plusOneHistoryMask_U;

                }
            }
            historyBuffer.process(errorBuffer.getWriteErrors(), bitWriter);
            errorBuffer.swap();
        }
    }

    private void decodeTail(@Unsigned int sets_U, @Unsigned byte[] soft_U) {
        // flush state registers
        // now we only shift in 0s, skipping 1-successors
        int highbit_U = 1 << order_U - 1;

        for (int i = (int) (Integer.toUnsignedLong(sets_U) - order_U + 1); Integer.compareUnsigned(i, sets_U) < 0; i++) {
            // lasterrors are the aggregate bit errors for the states of shiftregister for the previous
            // time slice
            @Unsigned byte[] history_U = historyBuffer.getSlice();
            // calculate the distance from all output states to our sliced bits
            if (soft_U != null) {
                if (softMeasurement == CORRECT_SOFT_LINEAR) {
                    for (int j = 0; Integer.compareUnsigned(j, 1 << rate_U) < 0; j++) {
                        distances_U[j] = Metric.softDistanceLinear(j, soft_U, rate_U, i * rate_U);
                    }
                } else {
                    for (int j = 0; Integer.compareUnsigned(j, 1 << rate_U) < 0; j++) {
                        distances_U[j] = (short) Metric.softDistanceQuadratic(j, soft_U, rate_U, i * rate_U);
                    }
                }
            } else {
                int out_U = bitReader.read(rate_U);
                for (int i2 = 0; Integer.compareUnsigned(i2, 1 << rate_U) < 0; i2++) {
                    distances_U[i2] = distance(i2, out_U);
                }
            }

            // a mask to get the high order bit from the shift register
            int numIter_U = highbit_U << 1;
            int skip_U = 1 << order_U - Integer.toUnsignedLong(sets_U - i);
            int baseSkip_U = skip_U >>> 1;

            int highbase_U = highbit_U >>> 1;
            for (int low_U = 0, high_U = highbit_U, base_U = 0; Integer.compareUnsigned(high_U, numIter_U) < 0; low_U += skip_U, high_U += skip_U, base_U += baseSkip_U) {
                int lowOutput_U = table_U[low_U];
                int highOutput_U = table_U[high_U];
                short lowDist_U = distances_U[lowOutput_U];
                short highDist_U = distances_U[highOutput_U];

                short lowPastError_U = errorBuffer.getReadError(base_U);
                short highPastError_U = errorBuffer.getReadError(highbase_U + base_U);

                short lowError_U = (short) (Short.toUnsignedInt(lowDist_U) + Short.toUnsignedInt(lowPastError_U));
                short highError_U = (short) (Short.toUnsignedInt(highDist_U) + Short.toUnsignedInt(highPastError_U));

                int successor_U = low_U;
                short error_U;
                byte historyMask_U;
                if (Short.toUnsignedInt(lowError_U) <= Short.toUnsignedInt(highError_U)) {
                    error_U = lowError_U;
                    historyMask_U = 0;
                } else {
                    error_U = highError_U;
                    historyMask_U = 1;
                }
                errorBuffer.setWriteError(successor_U, error_U);
                history_U[successor_U] = historyMask_U;

            }
            historyBuffer.processSkip(errorBuffer.getWriteErrors(), bitWriter, skip_U);
            errorBuffer.swap();
        }

    }

    public void fillTable(@Unsigned int rate_U, @Unsigned int order_U, @Unsigned short[] poly_U) {
        for(int i = 0; Integer.compareUnsigned(i, 1 << order_U) < 0; i++) {
            int out_U = 0;
            int mask_U = 1;
            for(int j_U = 0; Integer.compareUnsigned(j_U, rate_U) < 0; j_U++) {
                out_U |= Integer.bitCount(i & Short.toUnsignedInt(poly_U[j_U])) % 2 != 0 ? mask_U : 0;
                mask_U <<= 1;
            }
            table_U[i] = out_U;
        }
    }


}
