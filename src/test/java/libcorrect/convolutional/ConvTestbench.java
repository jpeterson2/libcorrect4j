/*
 * libcorrect4j
 * ConvTestbench.java
 * Created from include/correct/util/error-sim.h
 *              util/error-sim.c
 *              tests/convolutional.c               @ https://github.com/quiet/libcorrect
 */

package libcorrect.convolutional;

import java.util.Arrays;
import java.util.Random;

public class ConvTestbench {
    public final static int maxBlockLen = 4_096;
    public final static Random RANDOM = new Random(1);

    private byte[] encoded;
    private double[] v;
    private double[] corrupted;
    private byte[] soft;
    private double[] noise;
    private long enclen;
    private long enclenBytes;
    private long msgLen;

    private Convolutional conv;

    public boolean assertTestResult(Convolutional c, long testLength, long rate, long order, double ebN0, double errorRate) {
        boolean res = true;
        conv = c;
        msgLen = testLength;
        double bpskVoltage = 1.0 / Math.sqrt(2.0);
        double bpskSymEnergy = 2 * Math.pow(bpskVoltage, 2.0);
        double bpskBitEnergy = bpskSymEnergy * Double.parseDouble(Long.toUnsignedString(rate));

        long errorCount_U = testConv(ebN0, bpskBitEnergy, bpskVoltage);
        double observedErrorRate = Double.parseDouble(Long.toUnsignedString(errorCount_U)) / (Double.parseDouble(Long.toUnsignedString(testLength)) * 8);
        if(observedErrorRate > errorRate) {
            System.out.printf("Test failed, expected error rate=%.2e, observed error rate=%.2e @%.1fdB for rate %d order %d\n", errorRate, observedErrorRate, ebN0, rate, order);
            res = false;
        } else {
            System.out.printf("Test passed, expected error rate=%.2e, observed error rate=%.2e @%.1fdB for rate %d order %d\n", errorRate, observedErrorRate, ebN0, rate, order);
        }
        return res;
    }

    private long testConv(double ebN0, double bpskBitEnergy, double bpskVoltage) {
        long numErrors = 0;

        while(msgLen != 0) {
            long blockLen_U = Long.compareUnsigned(maxBlockLen, msgLen) < 0 ? maxBlockLen : msgLen;
            msgLen -= blockLen_U;

            byte[] msg_U = new byte[(int)blockLen_U];
            for(int j_U = 0; Long.compareUnsigned(Integer.toUnsignedLong(j_U), blockLen_U) < 0; j_U++) {
                msg_U[j_U] = (byte)(RANDOM.nextInt() % 256);
            }

            resize(blockLen_U);

            ErrorSim.buildWhiteNoise(noise, enclen, ebN0, bpskBitEnergy);
            numErrors += testConvNoise(msg_U, blockLen_U, bpskVoltage);
        }
        return numErrors;
    }

    public int testConvNoise(byte[] msg, long nBytes, double bpskVoltage) {
        enclen = conv.encodeLen(msg.length);
        encoded = conv.encode(msg);

        ErrorSim.encodeBpsk(encoded, v, enclen, bpskVoltage);

        corrupted = Arrays.copyOf(v, (int) enclen);
        ErrorSim.addWhiteNoise(corrupted, noise, enclen);
        ErrorSim.decodeBpskSoft(corrupted, soft, enclen, bpskVoltage);

        byte[] msgOut = conv.decodeSoft(soft, enclen);

        if(msgOut.length != msg.length) {
            System.out.printf("Expected to decode %d bytes, decoded %d bytes instead\n", msg.length, msgOut.length);
            return -1;
        }

        return (int) ErrorSim.distance(msg, msgOut, msg.length);
    }

    public void resize(long msgLen) {
        enclen = conv.encodeLen(msgLen);
        enclenBytes = Long.remainderUnsigned(enclen, 8) != 0 ?
                                        Long.divideUnsigned(enclen, 8) + 1 :
                                        Long.divideUnsigned(enclen, 8);
        if (encoded != null) {
            encoded = Arrays.copyOf(encoded, (int) enclenBytes);
        } else {
            encoded = new byte[(int) enclenBytes];
        }
        if (v != null) {
            v = Arrays.copyOf(v, (int) enclen);
        } else {
            v = new double[(int) enclen];
        }
        if (corrupted != null) {
            corrupted = Arrays.copyOf(corrupted, (int) enclen);
        } else {
            corrupted = new double[(int) enclen];
        }
        if (noise != null) {
            noise = Arrays.copyOf(noise, (int) enclen);
        } else {
            noise = new double[(int) enclen];
        }
        if (soft != null) {
            soft = Arrays.copyOf(soft, (int) enclen);
        } else {
            soft = new byte[(int) enclen];
        }
    }

}
