package org.pla.compression.util.Encoding;

public class DPValue {
    private double value;
    private int index;

    public DPValue(double value, int index) {
        this.value = value;
        this.index = index;
    }

    public double getValue() {
        return value;
    }

    public int getIndex() {
        return index;
    }
}