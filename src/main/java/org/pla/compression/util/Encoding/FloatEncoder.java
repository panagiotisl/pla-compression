package org.pla.compression.util.Encoding;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;

public class FloatEncoder {
    public static void write(float number, ByteArrayOutputStream outputStream) throws IOException {
        int intBits = Float.floatToIntBits(number);
        IntEncoder.write(intBits, outputStream);
    }

    public static float read(ByteArrayInputStream inputStream) throws IOException {
        int number = IntEncoder.read(inputStream);
        return Float.intBitsToFloat(number);
    }
}
