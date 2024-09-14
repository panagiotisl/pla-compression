package org.pla.compression.encodings;

public class Result {
    private byte[] binary;
    private int segmentsSize;

    public Result(byte[] byteArrayImproved, int size) {
        this.binary = byteArrayImproved;
        this.segmentsSize = size;
    }

    public byte[] getBinary() {
        return binary;
    }

    public int getSegmentsSize() {
        return segmentsSize;
    }
}
