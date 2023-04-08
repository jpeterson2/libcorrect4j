/*
 * libcorrect4j
 * Polynomial.java
 * Created from src/correct/reed-solomon/encode.c
 *              src/correct/reed-solomon/encode.c @ https://github.com/quiet/libcorrect
 */

package libcorrect.reed_solomon;

import java.util.Arrays;
import org.checkerframework.checker.signedness.qual.Unsigned;

import static libcorrect.reed_solomon.Polynomial.*;

public class ReedSolomon {
    public static final short correctRsPrimitivePolynomial_8_4_3_2_0 = 0x11d;        // x^8 + x^4 + x^3 + x^2 + 1
    public static final short correctRsPrimitivePolynomial_8_5_3_1_0 = 0x12b;        // x^8 + x^5 + x^3 + x + 1
    public static final short correctRsPrimitivePolynomial_8_5_3_2_0 = 0x12d;        // x^8 + x^5 + x^3 + x^2 + 1
    public static final short correctRsPrimitivePolynomial_8_6_3_2_0 = 0x14d;        // x^8 + x^6 + x^3 + x^2 + 1
    public static final short correctRsPrimitivePolynomial_8_6_4_3_2_1_0 = 0x15f;    // x^8 + x^6 + x^4 + x^3 + x^2 + x + 1;
    public static final short correctRsPrimitivePolynomial_8_6_5_1_0 = 0x163;        // x^8 + x^6 + x^5 + x + 1
    public static final short correctRsPrimitivePolynomial_8_6_5_2_0 = 0x165;        // x^8 + x^6 + x^5 + x^2 + 1
    public static final short correctRsPrimitivePolynomial_8_6_5_3_0 = 0x169;        // x^8 + x^6 + x^5 + x^3 + 1
    public static final short correctRsPrimitivePolynomial_8_6_5_4_0 = 0x171;        // x^8 + x^6 + x^5 + x^4 + 1
    public static final short correctRsPrimitivePolynomial_8_7_2_1_0 = 0x187;        // x^8 + x^7 + x^2 + x + 1
    public static final short correctRsPrimitivePolynomial_8_7_3_2_0 = 0x18d;        // x^8 + x^7 + x^3 + x^2 + 1
    public static final short correctRsPrimitivePolynomial_8_7_5_3_0 = 0x1a9;        // x^8 + x^7 + x^5 + x^3 + 1
    public static final short correctRsPrimitivePolynomial_8_7_6_1_0 = 0x1c3;        // x^8 + x^7 + x^6 + x + 1
    public static final short correctRsPrimitivePolynomial_8_7_6_3_2_1_0 = 0x1cf;    // x^8 + x^7 + x^6 + x^3 + x^2 + x + 1
    public static final short correctRsPrimitivePolynomial_8_7_6_5_2_1_0 = 0x1e7;    // x^8 + x^7 + x^6 + x^5 + x^2 + x + 1
    public static final short correctRsPrimitivePolynomial_8_7_6_5_4_2_0 = 0x1f5;    // x^8 + x^7 + x^6 + x^5 + x^4 + x^2 + 1
    public static final short correctRsPrimitivePolynomialCcsds      = 0x187;        // x^8 + x^7 + x^2 + x + 1



    private final long blockLength; // final assigned as constant - should it be glb or w/e? or is there something about subclassing that means we don't actually know enough about it's value
    private final long maxMessageLength;
    private final long minDistance;
    private final byte fConsecutiveRoot;
    private final @Unsigned byte gRootGap;
    private final Field field;
    private final Polynomial generator;
    private final byte[] generatorRoots;
    private byte[][] generatorRootExp;
    private final Polynomial encodedPolynomial;
    private final Polynomial encodedRemainder;
    private byte[] syndromes;
    private byte[] modifiedSyndromes;
    private Polynomial receivedPolynomial;
    private Polynomial errorLocator;
    private Polynomial errorLocatorLog;
    private Polynomial erasureLocator;
    private byte[] errorRoots;
    private byte[] errorVals;
    private byte[] errorLocations;

    private byte[][] elementExp;
    //  scratch (do no allocations at steady state)
    //  used during find_error_locator
    private Polynomial lastErrorLocator;

    // used during error value search
    private Polynomial errorEvaluator;
    private Polynomial errorLocatorDerivative;
    private final Polynomial[] initFromRootsScratch = new Polynomial[2];
    private boolean hasInitDecode;

    /**
     *  correctReedSolomon allocates and initializes an encoder/decoder for a given
     *  reed solomon error correction code. The block size must be 255 bytes with 8-bit symbols.
     * @param primitivePolynomial           Should be one of the given values in this file
     * @param firstConsecutiveRoot
     * @param generatorRootGap              Sane values for firstConsecutiveRoot and generatorRootGap
     *                                      are 1 and 1. Not all combinations of values produce valid codes.
     * @param numRoots                      Handle as many as numRoots/2 bytes having corruption
     *                                      and still recover the encoded payload. However, using
     *                                      more numRoots adds more parity overhead and substantially
     *                                      increases the computational time for decoding.
     */

    public ReedSolomon(short primitivePolynomial, byte firstConsecutiveRoot, @Unsigned byte generatorRootGap, long numRoots) {
        field = new Field(primitivePolynomial);
        blockLength = 255;
        minDistance = numRoots;

        maxMessageLength = blockLength - minDistance;
        fConsecutiveRoot = firstConsecutiveRoot;
        gRootGap = generatorRootGap;
        generatorRoots = new byte[(int) minDistance];

        generator = reedSolomonBuildGenerator((int) minDistance, generatorRoots);

        encodedPolynomial = new Polynomial((int) (blockLength - 1));
        encodedRemainder = new Polynomial((int) (blockLength - 1));

        hasInitDecode = false;
    }

    /**
     * encode uses the rs instance to encode parity information onto a block of data.
     * message length should be no more than the payload size for one block e.g. no more
     * (255 - min. distance). Shorter blocks will be encoded with virtual padding
     * where the padding is not emitted.
     * @param msg    a message to encode
     * @return       encoded message
     * @throws IllegalArgumentException  if message length is larger then (255 - min. distance)
     */
    public byte[] encode(byte[] msg) throws IllegalArgumentException {

        long msgLength = msg.length; // should be safe because length of array is nonnegative, we'll see if checker knows that :)
        byte[] encoded = new byte[(int) (msgLength+minDistance)];

        if (msgLength > maxMessageLength) {
            throw new IllegalArgumentException("ReedSolomon.encode: message length must be smaller than block length - min. distance");
        }

        long padLength = maxMessageLength - msgLength;

        // 0-fill the rest of the coefficients -- this length will always be > 0
        // because the order of this poly is block_length and the msg_length <= message_length
        // e.g. 255 and 223
        encodedPolynomial.flushCoeff();

        for (int i = 0; Long.compareUnsigned(Integer.toUnsignedLong(i), msgLength) < 0; i++) {
            // message goes from high order to low order but libcorrect polynomials go low to high
            // so we reverse on the way in and on the way out
            // we'd have to do a copy anyway so this reversal should be free
            encodedPolynomial.setCoeff((int) (encodedPolynomial.getOrder() - (i + padLength)), msg[i]); // I think Integer.toUnsignedLong should be able to consume @Unsigned
        }

        Polynomial.mod(field, encodedPolynomial, generator, encodedRemainder);

        // now return byte order to highest order to lowest order
        for (int i = 0; i < msgLength; i++) {
            encoded[i] = encodedPolynomial.getCoeff((int) (encodedPolynomial.getOrder() - (i + padLength))); // toUnsignedLong should be annotated such that arguments can be whatever signedness i think
        }

        for (int i = 0; i < minDistance; i++) {
            encoded[(int) (msgLength + Integer.toUnsignedLong(i))] = encodedRemainder.getCoeff((int) (minDistance - (i + 1)));
        }

        return encoded;
    }

    /**
     * decode uses the rs instance to decode a payload from a block containing
     * payload and parity bytes. This function can recover in spite of some bytes
     * being corrupted. In most cases, if the block is too corrupted, this function
     * will return -1 and not perform decoding. It is possible but unlikely that
     * the payload written to msg will contain errors when this function returns
     * a positive value.
     * @param encoded       encoded message
     * @return              decoded message or null if the block is too corrupted and cannot be recovered
     * @throws IllegalArgumentException  if encoded message length is larger than block length
     */
    public byte[] decode(byte[] encoded) throws IllegalArgumentException {

        long encodedLength = encoded.length;

        if (encodedLength > blockLength) {
            throw new IllegalArgumentException("ReedSolomon.decode: encoded message length must be smaller than block length");
        }

        // the message is the non-remainder part
        long msgLength = encodedLength - minDistance;
        byte[] msg = new byte[(int) msgLength];
        // if they handed us a nonfull block, we'll write in 0s
        long padLength = blockLength - encodedLength;

        if (!hasInitDecode) {
            // initialize rs for decoding
            createDecoder();
        }

        // we need to copy to our local buffer
        // the buffer we're given has the coordinates in the wrong direction
        // e.g. byte 0 corresponds to the 254th order coefficient
        // so we're going to flip and then write padding
        // the final copied buffer will look like
        // | rem (rs->min_distance) | msg (msg_length) | pad (pad_length) |

        for (int i = 0; i < encodedLength; i++) {
            receivedPolynomial.setCoeff(i, encoded[(int) (encodedLength - (i + 1))]);
        }

        // fill the pad_length with 0s
        for (int i = 0; i < padLength; i++) {
            receivedPolynomial.setCoeff((int) (Integer.toUnsignedLong(i) + encodedLength), (byte) 0);
        }


        boolean allZero = findSyndromes(receivedPolynomial);

        if (allZero) {
            // syndromes were all zero, so there was no error in the message
            // copy to msg and we are done
            for (int i = 0; i < msgLength; i++) {
                msg[i] = receivedPolynomial.getCoeff((int) (encodedLength - Integer.toUnsignedLong(i + 1)));
            }
            return msg;
        }
        int order = findErrorLocator(0);
        // XXX fix this vvvv
        errorLocator.setOrder(order);

        for (int i = 0; Integer.compare(i, errorLocator.getOrder()) <= 0; i++) {
            // this is a little strange since the coeffs are logs, not elements
            // also, we'll be storing log(0) = 0 for any 0 coeffs in the error locator
            // that would seem bad but we'll just be using this in chien search, and we'll skip all 0 coeffs
            // (you might point out that log(1) also = 0, which would seem to alias. however, that's ok,
            //   because log(1) = 255 as well, and in fact that's how it's represented in our log table)
            errorLocatorLog.setCoeff(i, field.log(Byte.toUnsignedInt(errorLocator.getCoeff(i))));
        }
        errorLocatorLog.setOrder(errorLocator.getOrder());
        if (!factorizeErrorLocator(0, errorLocatorLog, errorRoots, elementExp)) {
            // roots couldn't be found, so there were too many errors to deal with
            // RS has failed for this message
            return null;
        }

        findErrorLocations(errorLocator.getOrder());
        findErrorValues();

        for (int i = 0; Integer.compare(i, errorLocator.getOrder()) < 0; i++) {
            receivedPolynomial.setCoeff(Byte.toUnsignedInt(errorLocations[i]),
                    field.fieldSub(receivedPolynomial.getCoeff(Byte.toUnsignedInt(errorLocations[i])), errorVals[i]));
        }

        for (int i = 0; i < msgLength; i++) {
            msg[i] = receivedPolynomial.getCoeff((int) (encodedLength - Integer.toUnsignedLong(i + 1)));
        }
        return msg;
    }

     /**
     * decodeWithErasures decodes a payload from a block containing payload
     * and parity bytes. Additionally, the user can provide the indices of bytes
     * which have been suspected to be corrupted.
     * This erasure information is typically provided by a demodulating or receiving device.
     * This function can recover with some additional errors on top of the erasures.
     *
     * In order to successfully decode, the quantity (numErasures + 2*numErrors) must be less than
     * numRoots.
     *
     * erasure_locations should contain erasure_length items.
     * erasure_length should not exceed the number of parity
     * bytes encoded into this block.
     *
     * In most cases, if the block is too corrupted, this function
     * will return -1 and not perform decoding. It is possible but
     * unlikely that the payload written to msg will contain
     * errors when this function returns a positive value.
     *
     * msg should be long enough to contain a decoded payload for
     * this encoded block.
     * @return a positive number of bytes written to msg if it has decoded or -1 if it has encountered an error.
     */

    /**
     * decodeWithErasures decodes a payload from a block containing payload
     * and parity bytes. Additionally, the user can provide the indices of bytes
     * which have been suspected to be corrupted.
     * This erasure information is typically provided by a demodulating or receiving device.
     * This function can recover with some additional errors on top of the erasures.
     *
     * In order to successfully decode, the quantity (numErasures + 2*numErrors) must be less than
     * numRoots.
     * In most cases, if the block is too corrupted, this function will return -1 and not perform
     * decoding. It is possible but unlikely that the payload written to msg will contain
     * errors when this function returns a positive value.
     *
     * @param encoded                   encoded message
     * @param erasureLocations          erasure locations
     * @return                          decoded message or null if the block is too corrupted and cannot be recovered
     * @throws IllegalArgumentException if encoded message length is larger than block length or number of erasures is
     *                                  larger than min. distance
     */
    public byte[] decodeWithErasures(byte[] encoded, byte[] erasureLocations) throws IllegalArgumentException {
        long erasureLength = erasureLocations.length;
        if (erasureLength == 0) {
            return decode(encoded);
        }

        long encodedLength = encoded.length;
        if (encodedLength > blockLength) {
            throw new IllegalArgumentException("ReedSolomon.decodeWithErasures: encoded message length must be smaller than block length");
        }
        if (erasureLength > minDistance) {
            throw new IllegalArgumentException("ReedSolomon.decodeWithErasures: erasures  length must be smaller than min distance");
        }

        // the message is the non-remainder part
        long msgLength = encodedLength - minDistance;
        // if they handed us a nonfull block, we'll write in 0s
        long padLength = blockLength - encodedLength;

        byte[] msg = new byte[(int) msgLength];

        if (!hasInitDecode) {
            createDecoder();
        }

        // we need to copy to our local buffer
        // the buffer we're given has the coordinates in the wrong direction
        // e.g. byte 0 corresponds to the 254th order coefficient
        // so we're going to flip and then write padding
        // the final copied buffer will look like
        // | rem (rs->min_distance) | msg (msg_length) | pad (pad_length) |

        for (int i = 0; i < encodedLength; i++) {
            receivedPolynomial.setCoeff(i, encoded[(int) (encodedLength - Integer.toUnsignedLong(i + 1))]);
        }

        // fill the pad_length with 0s
        for (int i = 0; i < padLength; i++) {
            receivedPolynomial.setCoeff((int) (i + encodedLength), (byte) 0);
        }

        for (int i = 0; i < erasureLength; i++) {
            // remap the coordinates of the erasures
            errorLocations[i] = (byte) (blockLength - (Byte.toUnsignedLong(erasureLocations[i]) + padLength + 1));
        }

        findErrorRootsFromLocations(gRootGap, (int) erasureLength);

        findErrorLocatorFromRoots((int) erasureLength, erasureLocator, initFromRootsScratch);

        boolean allZero = findSyndromes(receivedPolynomial);

        if (allZero) {
            // syndromes were all zero, so there was no error in the message
            // copy to msg and we are done
            for (int i = 0; i < msgLength; i++) {
                msg[i] = receivedPolynomial.getCoeff((int) (encodedLength - Integer.toUnsignedLong(i + 1)));
            }
            return msg;
        }

        findModifiedSyndromes(erasureLocator);

        byte[] syndromeCopy_U = Arrays.copyOf(syndromes, (int) minDistance);

        for (int i = (int) erasureLength; i < minDistance; i++) {
            syndromes[(int) (Integer.toUnsignedLong(i) - erasureLength)] = modifiedSyndromes[i];
        }

        int order = findErrorLocator(erasureLength);
        // XXX fix this vvvv
        errorLocator.setOrder(order);

        for (int i = 0; i <= errorLocator.getOrder(); i++) {
            // this is a little strange since the coeffs are logs, not elements
            // also, we'll be storing log(0) = 0 for any 0 coeffs in the error locator
            // that would seem bad but we'll just be using this in chien search, and we'll skip all 0 coeffs
            // (you might point out that log(1) also = 0, which would seem to alias. however, that's ok,
            //   because log(1) = 255 as well, and in fact that's how it's represented in our log table)
            errorLocatorLog.setCoeff(i, field.log((Byte.toUnsignedInt(errorLocator.getCoeff(i)))));
        }
        errorLocatorLog.setOrder(errorLocator.getOrder());

        if (!factorizeErrorLocator((int) erasureLength, errorLocatorLog, errorRoots, elementExp)) {
            // roots couldn't be found, so there were too many errors to deal with
            // RS has failed for this message
            return null;
        }

        Polynomial tempPoly = new Polynomial((int) (Integer.toUnsignedLong(errorLocator.getOrder()) + erasureLength)); // I believe sign extension is an issue here -  therefore the annotation on toUnsignedLong is incorrect i think
        mul(field, erasureLocator, errorLocator, tempPoly);
        Polynomial placeholderPoly = errorLocator;
        errorLocator = tempPoly;

        findErrorLocations(errorLocator.getOrder());

        syndromes = Arrays.copyOf(syndromeCopy_U, (int) minDistance);
        findErrorValues();

        for (int i = 0; i < errorLocator.getOrder(); i++) {
            receivedPolynomial.setCoeff(Byte.toUnsignedInt(errorLocations[i]),
                    field.fieldSub(receivedPolynomial.getCoeff(Byte.toUnsignedInt(errorLocations[i])), errorVals[i]));
        }

        errorLocator = placeholderPoly;
        for (int i = 0; i < msgLength; i++) {
            msg[i] = receivedPolynomial.getCoeff((int) (encodedLength - Integer.toUnsignedLong(i + 1)));
        }
        return msg;
    }

    /**
     * Calculate all syndromes of the received polynomial at the roots of the generator
     * because we're evaluating at the roots of the generator, and because the transmitted
     * polynomial was made to be a product of the generator, we know that the transmitted
     * polynomial is 0 at these roots
     * any nonzero syndromes we find here are the values of the error polynomial evaluated
     * at these roots, so these values give us a window into the error polynomial. if
     * these syndromes are all zero, then we can conclude the error polynomial is also
     * zero. if they're nonzero, then we know our message received an error in transit.
     * @param msgpoly       Received polynomial
     * @return              True if syndromes are all zero, false otherwise
     */
    private boolean findSyndromes(Polynomial msgpoly) {
        boolean allZero = true;
        Arrays.fill(syndromes, 0, (int) minDistance, (byte) 0);
        for (int i = 0; i < minDistance; i++) {
            // profiling reveals that this function takes about 50% of the cpu time of
            // decoding. so, in order to speed it up a little, we precompute and save
            // the successive powers of the roots of the generator, which are
            // located in generator_root_exp
            byte eval_U = msgpoly.evalLut(field, generatorRootExp[i]);
            if (eval_U != 0) {
                allZero = false;
            }
            syndromes[i] = eval_U;
        }
        return allZero;
    }

    /**
     * Berlekamp-Massey algorithm to find LFSR that describes syndromes
     * @param numErasures
     * @return Returns number of errors and writes the error locator polynomial to errorLocator
     */
    private int findErrorLocator(long numErasures) {
        int numerrors = 0;
        errorLocator.flushCoeff();

        // initialize to f(x) = 1
        errorLocator.setCoeff(0, (byte) 1);
        errorLocator.setOrder(0);
        lastErrorLocator.setOrder(0);

        lastErrorLocator.copyCoeff(errorLocator);

        byte discrepancy;
        byte lastDiscrepancy = 1;
        int delayLength = 1;

        for (int i = errorLocator.getOrder(); i < (minDistance - numErasures); i++) {
            discrepancy = syndromes[i];
            for (int j = 1; j <= numerrors; j++) {
                discrepancy = (byte) field.fieldAdd(discrepancy, field.fieldMul(errorLocator.getCoeff(j), syndromes[i - j]));
            }

            if (discrepancy == 0) {
                // our existing LFSR describes the new syndrome as well
                // leave it as-is but update the number of delay elements
                //   so that if a discrepancy occurs later we can eliminate it
                delayLength++;
                continue;
            }

            if (2 * numerrors <= i) {
                // there's a discrepancy, but we still have room for more taps
                // lengthen LFSR by one tap and set weight to eliminate discrepancy

                // shift the last locator by the delay length, multiply by discrepancy,
                //   and divide by the last discrepancy
                // we move down because we're shifting up, and this prevents overwriting
                for (int j = lastErrorLocator.getOrder(); j >= 0; j--) { // TODO: rewrite with compareUnsigned as intended
                    // the bounds here will be ok since we have a headroom of numerrors
                    lastErrorLocator.setCoeff(j + delayLength, field.fieldDiv(
                            field.fieldMul(lastErrorLocator.getCoeff(j), discrepancy), lastDiscrepancy));
                }
                for (int j = delayLength - 1; j >= 0; j--) {
                    lastErrorLocator.setCoeff(j, (byte) 0);
                }

                // locator = locator - last_locator
                // we will also update last_locator to be locator before this loop takes place
                byte temp;
                for (int j = 0; j <= (lastErrorLocator.getOrder() + delayLength); j++) {
                    temp = errorLocator.getCoeff(j);
                    errorLocator.setCoeff(j, field.fieldAdd(errorLocator.getCoeff(j), lastErrorLocator.getCoeff(j)));
                    lastErrorLocator.setCoeff(j, temp);
                }
                int tempOrder = errorLocator.getOrder();
                errorLocator.setOrder(lastErrorLocator.getOrder() + delayLength);
                lastErrorLocator.setOrder(tempOrder);

                // now last_locator is locator before we started,
                //   and locator is (locator - (discrepancy/last_discrepancy) * x^(delay_length) * last_locator)

                numerrors = i + 1 - numerrors;
                lastDiscrepancy = discrepancy;
                delayLength = 1;
                continue;
            }

            // no more taps
            // unlike the previous case, we are preserving last locator,
            //    but we'll update locator as before
            // we're basically flattening the two loops from the previous case because
            //    we no longer need to update last_locator
            for (int j = lastErrorLocator.getOrder(); j >= 0; j--) {
                errorLocator.setCoeff(j + delayLength, (byte) field.fieldAdd(errorLocator.getCoeff(j + delayLength),
                        field.fieldDiv(field.fieldMul(lastErrorLocator.getCoeff(j), discrepancy), lastDiscrepancy)));
            }
            errorLocator.setOrder((lastErrorLocator.getOrder() + delayLength) > errorLocator.getOrder() ?
                    lastErrorLocator.getOrder() + delayLength : errorLocator.getOrder());
            delayLength++;
        }
        return errorLocator.getOrder();
    }

    /**
     * Find the roots of the error locator polynomial (Chien search)
     * @param numSkip
     * @param locatorLog
     * @param roots
     * @param elementExp
     * @return  True if errors are recoverable, false otherwise
     */
    public boolean factorizeErrorLocator(int numSkip, Polynomial locatorLog, byte[] roots, byte[][] elementExp) {
        // normally it'd be tricky to find all the roots
        // but, the finite field is awfully finite...
        // just brute force search across every field element
        int root = numSkip;
        for (int i = 0; i < locatorLog.getOrder(); i++) {
            roots[numSkip + i] = 0;
        }
        for (short i = 0; Short.toUnsignedInt(i) < 256; i++) {
            // we make two optimizations here to help this search go faster
            // a) we have precomputed the first successive powers of every single element
            //   in the field. we need at most n powers, where n is the largest possible
            //   degree of the error locator
            // b) we have precomputed the error locator polynomial in log form, which
            //   helps reduce some lookups that would be done here
            if (locatorLog.evalLogLut(field, elementExp[Short.toUnsignedInt(i)]) == 0) {
                roots[root] = (byte) i;
                root++;
            }
        }
        // this is where we find out if we are have too many errors to recover from
        // berlekamp-massey may have built an error locator that has 0 discrepancy
        // on the syndromes but doesn't have enough roots
        return root == locatorLog.getOrder() + numSkip;
    }

    /**
     * use error locator and syndromes to find the error evaluator polynomial
     * @param locator
     * @param syndromes
     * @param errorEvaluator
     */
    public void findErrorEvaluator(Polynomial locator, Polynomial syndromes, Polynomial errorEvaluator) {
        // the error evaluator, omega(x), is S(x)*Lamba(x) mod x^(2t)
        // where S(x) is a polynomial constructed from the syndromes
        //   S(1) + S(2)*x + ... + S(2t)*x(2t - 1)
        // and Lambda(x) is the error locator
        // the modulo is implicit here -- we have limited the max length of error_evaluator,
        //   which polynomial_mul will interpret to mean that it should not compute
        //   powers larger than that, which is the same as performing mod x^(2t)
        Polynomial.mul(field, locator, syndromes, errorEvaluator);
    }


    /**
     * use error locator, error roots and syndromes to find the error values
     * that is, the elements in the finite field which can be added to the received
     *   polynomial at the locations of the error roots in order to produce the
     *   transmitted polynomial
     * forney algorithm
     */
    public void findErrorValues() {
        // error value e(j) = -(X(j)^(1-c) * omega(X(j)^-1))/(lambda'(X(j)^-1))
        // where X(j)^-1 is a root of the error locator, omega(X) is the error evaluator,
        //   lambda'(X) is the first formal derivative of the error locator,
        //   and c is the first consecutive root of the generator used in encoding

        // first find omega(X), the error evaluator
        // we generate S(x), the polynomial constructed from the roots of the syndromes
        // this is *not* the polynomial constructed by expanding the products of roots
        // S(x) = S(1) + S(2)*x + ... + S(2t)*x(2t - 1)
        Polynomial syndromePoly = new Polynomial((int) (minDistance - 1), syndromes);
        errorEvaluator.flushCoeff();
        findErrorEvaluator(errorLocator, syndromePoly, errorEvaluator);

        // now find lambda'(X)
        errorLocatorDerivative.setOrder(errorLocator.getOrder() - 1);
        formalDerivative(field, errorLocator, errorLocatorDerivative);

        // calculate each e(j)
        for (int i = 0; i < errorLocator.getOrder(); i++) {
            if (Byte.toUnsignedInt(errorRoots[i]) == 0) {
                continue;
            }
            errorVals[i] = field.fieldMul(field.fieldPow(errorRoots[i], Byte.toUnsignedInt(fConsecutiveRoot) - 1),
                    field.fieldDiv(
                            errorEvaluator.evalLut(field, elementExp[Byte.toUnsignedInt(errorRoots[i])]),
                            errorLocatorDerivative.evalLut(field, elementExp[Byte.toUnsignedInt(errorRoots[i])])));
        }
    }

    /**
     *
     * @param numErrors
     */
    public void findErrorLocations(int numErrors) {
        for (int i = 0; i < numErrors; i++) {
            // the error roots are the reciprocals of the error locations, so div 1 by them

            // we do mod 255 here because the log table aliases at index 1
            // the log of 1 is both 0 and 255 (alpha^255 = alpha^0 = 1)
            // for most uses it makes sense to have log(1) = 255, but in this case
            // we're interested in a byte index, and the 255th index is not even valid
            // just wrap it back to 0

            if (Byte.toUnsignedInt(errorRoots[i]) == 0) {
                continue;
            }

            byte loc_U = field.fieldDiv((byte) 1, errorRoots[i]);
            for (int j = 0; Integer.compareUnsigned(j, 256) < 0; j++) {
                if (field.fieldPow((byte) j, gRootGap) == loc_U) {
                    errorLocations[i] = field.log(j);
                    break;
                }
            }
        }
    }

    /**
     * erasure method -- take given locations and convert to roots
     * this is the inverse of findErrorLocations
     * @param generatorRootGap
     * @param numErrors
     */
    private void findErrorRootsFromLocations(@Unsigned byte generatorRootGap, int numErrors) {
        for (int i = 0; i < numErrors; i++) {
            byte loc_U = field.fieldPow(field.exp(Byte.toUnsignedInt(errorLocations[i])), generatorRootGap);
            // field_element_t loc = field.exp[error_locations[i]];
            errorRoots[i] = field.fieldDiv((byte) 1, loc_U);
        }

    }

    /**
     * erasure method -- given the roots of the error locator, create the polynomial
     * @param numErrors
     * @param errorLocator
     * @param scratch
     * @return
     */
    private void findErrorLocatorFromRoots(int numErrors, Polynomial errorLocator, Polynomial[] scratch) {
        // multiply out roots to build the error locator polynomial
        errorLocator.initFromRoots(field, numErrors, errorRoots, scratch);
    }

    /**
     * erasure method
     * @param errorLocator
     */
    private void findModifiedSyndromes(Polynomial errorLocator) {
        Polynomial syndromePoly = new Polynomial((int) (minDistance - 1), syndromes);
        Polynomial modifiedSyndromePoly = new Polynomial((int) (minDistance - 1));

        mul(field, errorLocator, syndromePoly, modifiedSyndromePoly);
        modifiedSyndromes = Arrays.copyOf(modifiedSyndromePoly.getCoeff(),(int) minDistance);
    }

    /**
     * Creates decoder
     */
    private void createDecoder() {
        hasInitDecode = true;
        syndromes = new byte[(int) minDistance];
        modifiedSyndromes = new byte[(int) (2 * minDistance)];
        receivedPolynomial = new Polynomial((int) (blockLength - 1));
        errorLocator = new Polynomial((int) minDistance);
        errorLocatorLog = new Polynomial((int) minDistance);
        erasureLocator = new Polynomial((int) minDistance);


        errorRoots = new byte[(int) (2 * minDistance)];
        errorVals = new byte[(int) minDistance];
        errorLocations = new byte[(int) minDistance];

        lastErrorLocator = new Polynomial((int) minDistance);
        errorEvaluator = new Polynomial((int) minDistance -1);
        errorLocatorDerivative = new Polynomial((int) (minDistance - 1));

        // calculate and store the first block_length powers of every generator root
        // we would have to do this work in order to calculate the syndromes
        // if we save it, we can prevent the need to recalculate it on subsequent calls
        // total memory usage is min_distance * block_length bytes e.g. 32 * 255 ~= 8k
        generatorRootExp = new byte[(int) minDistance][];
        for (int i = 0; i < minDistance; i++) {
            generatorRootExp[i] = new byte[(int) blockLength];
            buildExpLut(field, generatorRoots[i], (int) (blockLength - 1), generatorRootExp[i]);
        }

        // calculate and store the first min_distance powers of every element in the field
        // we would have to do this for chien search anyway, and its size is only 256 * min_distance bytes
        // for min_distance = 32 this is 8k of memory, a pittance for the speedup we receive in exchange
        // we also get to reuse this work during error value calculation
        elementExp = new byte[256][];
        for (int i = 0; i < 256; i++) {
            elementExp[i] = new byte[(int) minDistance];
            buildExpLut(field, (byte) i, (int) (minDistance - 1), elementExp[i]);
        }

        initFromRootsScratch[0] = new Polynomial((int) minDistance);
        initFromRootsScratch[1] = new Polynomial((int) minDistance);

    }

    /**
     * Build generator from roots
     * coeff must be of size nroots + 1
     * e.g. 2 roots (x + alpha)(x + alpha^2) yields a poly with 3 terms x^2 + g0*x + g1
     * @param nroots
     * @param roots
     * @return
     */
    private Polynomial reedSolomonBuildGenerator(int nroots,  byte[] roots) {
        for (int i = 0; i < nroots; i++) {
            roots[i] = field.exp(Integer.remainderUnsigned(
                    gRootGap * (i + Byte.toUnsignedInt(fConsecutiveRoot)), 255));
        }
        return new Polynomial(field, nroots, roots);

    }

    /**
     * Debug print
     */
    public void debugPrint() {
        for (int i = 0; Integer.compareUnsigned(i, 256) < 0; i++) {
            System.out.printf("%3d  %3d    %3d  %3d\n", i, Byte.toUnsignedInt(field.exp(i)), i, Byte.toUnsignedInt(field.log(i)));
        }
        System.out.println();

        System.out.print("roots: ");
        for (int i = 0; i < minDistance; i++) {
            System.out.print(Byte.toUnsignedInt(generatorRoots[i]));
            if (i < (minDistance - 1)) {
                System.out.print(", ");
            }
        }
        System.out.println('\n');

        System.out.print("generator: ");
        for (int i = 0; i < (generator.getOrder() + 1); i++) {
            System.out.print(Byte.toUnsignedInt(generator.getCoeff(i)) + "*x^" + i);
            if (i < generator.getOrder()) {
                System.out.print(" + ");
            }
        }
        System.out.println('\n');

        System.out.print("generator (alpha format): ");
        for (int i = generator.getOrder() + 1; i > 0; i--) {
            System.out.print("alpha^" + Byte.toUnsignedInt(
                    field.log(Byte.toUnsignedInt(generator.getCoeff(i - 1)))) + "*x^" + (i - 1));
            if (i > 1) {
                System.out.print(" + ");
            }
        }
        System.out.println("\n");

        System.out.print("remainder: ");
        boolean hasPrinted = false;
        for (int i = 0; i < (encodedRemainder.getOrder() + 1); i++) {
            if (encodedRemainder.getCoeff(i) == 0) {
                continue;
            }
            if (hasPrinted) {
                System.out.print(" + ");
            }
            hasPrinted = true;
            System.out.print(Byte.toUnsignedInt(errorLocator.getCoeff(i)) + "*x^" + i);
        }
        System.out.println("\n");

        System.out.print("syndromes: ");
        for (int i = 0; i < minDistance; i++) {
            System.out.print(Byte.toUnsignedInt(syndromes[i]));
            if (i < (minDistance - 1)) {
                System.out.print(", ");
            }
        }
        System.out.println("\n");

        System.out.println("numerrors: " + errorLocator.getOrder() + "\n");

        System.out.print("error locator: ");
        hasPrinted = false;
        for (int i = 0; i < (errorLocator.getOrder() + 1); i++) {
            if (errorLocator.getCoeff(i) == 0) {
                continue;
            }
            if (hasPrinted) {
                System.out.print(" + ");
            }
            hasPrinted = true;
            System.out.print(Byte.toUnsignedInt(errorLocator.getCoeff(i)) + "*x^" + i);
        }
        System.out.println("\n");

        System.out.print("error roots: ");
        for (int i = 0; i < (errorLocator.getOrder()); i++) {
            System.out.print(eval(field, errorLocator, errorRoots[i]) + "@" + Byte.toUnsignedInt(errorRoots[i]));
            if (i < (errorLocator.getOrder() - 1)) {
                System.out.print(", ");
            }
        }
        System.out.println("\n");

        hasPrinted = false;
        for(int i = 0; i < errorEvaluator.getOrder(); i++) {
            if(errorEvaluator.getCoeff(i) == 0) {
                continue;
            }
            if(hasPrinted) {
                System.out.print(" + ");
            }
            hasPrinted = true;
            System.out.print(Byte.toUnsignedInt(errorEvaluator.getCoeff(i)) + "*x^" + i);
        }
        System.out.println("\n");

        hasPrinted = false;
        for(int i = 0; i < errorLocatorDerivative.getOrder(); i++) {
            if(errorLocatorDerivative.getCoeff(i) == 0) {
                continue;
            }
            if(hasPrinted) {
                System.out.print(" + ");
            }
            hasPrinted = true;
            System.out.print(Byte.toUnsignedInt(errorLocatorDerivative.getCoeff(i)) + "*x^" + i);
        }
        System.out.println("\n");

        System.out.print("error locator: ");
        for(int i = 0; i < errorLocator.getOrder(); i++) {
            System.out.print(Byte.toUnsignedInt(errorVals[i]) + "@" + Byte.toUnsignedInt(errorLocations[i]));
            if(i < (errorLocator.getOrder() - 1)) {
                System.out.print(", ");
            }
        }
        System.out.println("\n");

    }


}
