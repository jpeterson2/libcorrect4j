/*
 * libcorrect4j
 * ErrorBuffer.java
 * Created from src/correct/convolutional/error_buffer.c @ https://github.com/quiet/libcorrect
 */

package libcorrect.convolutional;

import java.util.Arrays;

import org.checkerframework.checker.signedness.qual.Unsigned;

public class ErrorBuffer {
    private int index;
    private final @Unsigned short[][] errors_U;
    private @Unsigned int numStates_U;

    public ErrorBuffer(@Unsigned int numStates_U) {
        // how large are the error buffers?
        this.numStates_U = numStates_U;

        // save two error metrics, one for last round and one for this
        // (double buffer)
        // the error metric is the aggregated number of bit errors found
        //   at a given path which terminates at a particular shift register state
        errors_U = new @Unsigned short[2][];
        this.errors_U[0] = new short[numStates_U];
        this.errors_U[1] = new short[numStates_U];

        // which buffer are we using, 0 or 1?
        this.index = 0;
    }

    public void reset() {
        Arrays.fill(errors_U[0], (short)0);
        Arrays.fill(errors_U[1], (short)0);
        index = 0;
    }
    public void swap() {
        index = (index+1)%2;
    }
    public @Unsigned short getReadError(int i) {
        return errors_U[index][i];
    }
    public void setReadError(int i, @Unsigned short v) {
        errors_U[index][i] = v;
    }
    public @Unsigned short getWriteError(int i) {
        return errors_U[(index+1)%2][i];
    }
    public void setWriteError(int i, @Unsigned short v) {
        errors_U[(index+1)%2][i] = v;
    }
    public @Unsigned short[] getWriteErrors() {
        return errors_U[(index+1)%2];
    }


}

