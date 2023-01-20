/*
 * libcorrect4j
 * BitReader.java
 * Created from src/correct/convolutional/bit.c @ https://github.com/quiet/libcorrect
 */

package libcorrect.convolutional;

import org.checkerframework.checker.signedness.qual.Unsigned;

public class BitReader {

	private static final @Unsigned byte[] reverseTable_U = createReverseTable();

	private static @Unsigned byte[] createReverseTable() {
        @Unsigned byte[] rT_U = new @Unsigned byte[256];
		for(int i = 0; i < 256; i++) {
			rT_U[i] =
        (byte)((i & 0x80) >> 7 | (i & 0x40) >> 5 | (i & 0x20) >> 3 | 
               (i & 0x10) >> 1 | (i & 0x08) << 1 | (i & 0x04) << 3 | 
               (i & 0x02) << 5 | (i & 0x01) << 7);
		}
		return rT_U;
	}

	@Unsigned byte currentByte_U;
	@Unsigned long byteIndex_U;
  	@Unsigned long len_U;
  	@Unsigned long currentByteLen_U;
 	@Unsigned byte[] bytes_U;

	public BitReader(@Unsigned byte[] bytes_U, @Unsigned long len_U) {
		if(bytes_U != null) {
			reconfigure(bytes_U, len_U);
		}
	}

	public void reconfigure(@Unsigned byte[] bytes_U, @Unsigned long len_U) {
		this.bytes_U = bytes_U;
		this.len_U = len_U;
		this.currentByteLen_U = 8;
		this.currentByte_U=bytes_U[0];
		this.byteIndex_U=0;
	}

	public @Unsigned byte read(@Unsigned int n_U) {
		int read_U = 0;
		int nCopy_U = n_U;

		if(Long.compareUnsigned(this.currentByteLen_U, Integer.toUnsignedLong(n_U)) < 0) {
			read_U = Byte.toUnsignedInt(this.currentByte_U) & (1 << this.currentByteLen_U) - 1;
			this.byteIndex_U++;
			this.currentByte_U = this.bytes_U[(int)this.byteIndex_U];
			n_U = (int)(Integer.toUnsignedLong(n_U) - this.currentByteLen_U);
			this.currentByteLen_U=8;
			read_U <<= n_U;
		}

		byte copyMask_U = (byte)((1 << n_U) - 1);
		copyMask_U = (byte)(Byte.toUnsignedInt(copyMask_U) << this.currentByteLen_U - Integer.toUnsignedLong(n_U));
		read_U |= (Byte.toUnsignedInt(this.currentByte_U) & Byte.toUnsignedInt(copyMask_U)) >> this.currentByteLen_U - Integer.toUnsignedLong(n_U);
		this.currentByteLen_U = this.currentByteLen_U - Integer.toUnsignedLong(n_U);
		return (byte)(Byte.toUnsignedInt(reverseTable_U[read_U]) >>> 8 - nCopy_U);
	}
}
