package org.pla.compression.encodings;

import com.github.luben.zstd.Zstd;
import me.lemire.integercompression.IntWrapper;
import me.lemire.integercompression.IntegerCODEC;
import me.lemire.integercompression.Simple16;
import org.pla.compression.encodings.encoders.FloatEncoder;
import org.pla.compression.encodings.encoders.UIntEncoder;
import org.pla.compression.encodings.encoders.VariableByteEncoder;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

/**
 * Mix-Piece Algorithm for Compressing Time-Series Data
 */
public class MixPiece {
    private static ArrayList<MixPieceSegment> perBSegments;
    private static ArrayList<MixPieceSegment> perASegments;
    private static ArrayList<MixPieceSegment> restSegments;

    private static double epsilon;
    private static int globalMinB;
    private static long lastTimeStamp;

//    private static IntegerCODEC CODEC = new Composition(new Simple16(), new VariableByte());
    private static IntegerCODEC CODEC = new Simple16();

    /**
     * Compress a list of Points and return a binary representation
     *
     * @param points       Time-series data
     * @param error        Maximum absolute error
     * @param variableByte
     * @param zstd
     * @return Binary representation
     * @throws Exception
     */
    public static Result compress(List<Point> points, double error, boolean variableByte, boolean zstd) throws Exception {
        if (points.isEmpty() || error <= 0) throw new Exception();

        epsilon = error;
        lastTimeStamp = points.get(points.size() - 1).getTimestamp();
        ArrayList<MixPieceSegment> segments = compress(points);
        merge(segments);
        return new Result(toByteArray(variableByte, zstd), segments.size());
    }

    public static Result compressTunablePeekAhead(List<Point> points, double error, boolean variableByte, boolean zstd, double pow) throws Exception {
        if (points.isEmpty() || error <= 0) throw new Exception();

        epsilon = error;
        lastTimeStamp = points.get(points.size() - 1).getTimestamp();
        ArrayList<MixPieceSegment> segments = compressTunablePeekAhead(points, pow);
        merge(segments);
        return new Result(toByteArrayImproved(variableByte, zstd), segments.size());
    }

    public static Result compressQuantOptimal(List<Point> points, double error, boolean variableByte, boolean zstd, double pow) throws Exception {
        if (points.isEmpty() || error <= 0) throw new Exception();

        epsilon = error;
        lastTimeStamp = points.get(points.size() - 1).getTimestamp();
        ArrayList<MixPieceSegment> segments = compressQuantOptimal(points, pow);
        merge(segments);
        return new Result(toByteArrayImproved(variableByte, zstd), segments.size());
    }

    /**
     * Decompress a binary representation and return a list of Points
     *
     * @param binary       Binary representation
     * @param variableByte
     * @param zstd
     * @return Time-series data
     */
    public static List<Point> decompress(byte[] binary, boolean variableByte, boolean zstd) {
        readByteArray(binary, variableByte, zstd);
        return toPoints();
    }

    public static List<Point> decompressImproved(byte[] binary, boolean variableByte, boolean zstd) {
        readByteArrayImproved(binary, variableByte, zstd);
        return toPoints();
    }

    private static double quantization(double value, int mode) {
        if (mode == 1) return (int) Math.ceil(value / epsilon) * epsilon;
        else if (mode == 2) return (int) Math.floor(value / epsilon) * epsilon;
        else return Math.round(value / epsilon) * epsilon;
    }

    private static int createSegment(int startIdx, List<Point> points, ArrayList<MixPieceSegment> segments, int quantizationMode) {
        long initTimestamp = points.get(startIdx).getTimestamp();
        double b = quantization(points.get(startIdx).getValue(), quantizationMode);
        if (startIdx + 1 == points.size()) {
            segments.add(new MixPieceSegment(initTimestamp, -Double.MAX_VALUE, Double.MAX_VALUE, b));
            return startIdx + 1;
        }
        double aMax = ((points.get(startIdx + 1).getValue() + epsilon) - b) / (points.get(startIdx + 1).getTimestamp() - initTimestamp);
        double aMin = ((points.get(startIdx + 1).getValue() - epsilon) - b) / (points.get(startIdx + 1).getTimestamp() - initTimestamp);
        if (startIdx + 2 == points.size()) {
            segments.add(new MixPieceSegment(initTimestamp, aMin, aMax, b));
            return startIdx + 2;
        }

        for (int idx = startIdx + 2; idx < points.size(); idx++) {
            double upValue = points.get(idx).getValue() + epsilon;
            double downValue = points.get(idx).getValue() - epsilon;

            double upLim = aMax * (points.get(idx).getTimestamp() - initTimestamp) + b;
            double downLim = aMin * (points.get(idx).getTimestamp() - initTimestamp) + b;
            if ((downValue > upLim || upValue < downLim)) {
                segments.add(new MixPieceSegment(initTimestamp, aMin, aMax, b));
                return idx;
            }

            if (upValue < upLim)
                aMax = Math.max((upValue - b) / (points.get(idx).getTimestamp() - initTimestamp), aMin);
            if (downValue > downLim)
                aMin = Math.min((downValue - b) / (points.get(idx).getTimestamp() - initTimestamp), aMax);
        }
        segments.add(new MixPieceSegment(initTimestamp, aMin, aMax, b));

        return points.size();
    }

    private static ArrayList<MixPieceSegment> compress(List<Point> points) {
        ArrayList<MixPieceSegment> segments = new ArrayList<>();
        int currentIdx = 0;
        while (currentIdx < points.size()) {
            int currentCeilIdx = createSegment(currentIdx, points, segments, 1);
            int currentFloorIdx = createSegment(currentIdx, points, segments, 2);
            if (currentCeilIdx > currentFloorIdx) {
                segments.remove(segments.size() - 1);
                currentIdx = currentCeilIdx;
            } else if (currentCeilIdx < currentFloorIdx) {
                segments.remove(segments.size() - 2);
                currentIdx = currentFloorIdx;
            } else {
                double firstValue = points.get(currentIdx).getValue();
                if (Math.round(firstValue / epsilon) == Math.ceil(firstValue / epsilon))
                    segments.remove(segments.size() - 1);
                else segments.remove(segments.size() - 2);
                currentIdx = currentFloorIdx;
            }
            globalMinB = (int) Math.min(globalMinB, segments.get(segments.size() - 1).getB() / epsilon);
        }

        return segments;
    }

    private static ArrayList<MixPieceSegment> compressTunablePeekAhead(List<Point> points, double pow) {
        ArrayList<MixPieceSegment> segments = new ArrayList<>();
        int startIdx = 0;
        while (startIdx < points.size()) {
            startIdx = Encoding.addSegment(startIdx, pow, points, epsilon, segments);
        }
        return segments;
    }

    private static ArrayList<MixPieceSegment> compressQuantOptimal(List<Point> points, double pow) {
        ArrayList<MixPieceSegment> segments = new ArrayList<>();
        Map<Integer, List<MixPieceSegment>> possibleMixPieceSegments = new TreeMap<>();
        for (int i = 0; i <= points.size() - 1; i++) {
            List<MixPieceSegment> segmentsFromStartIdx = Encoding.createMixPieceSegmentsFromStartIdx(i, points, epsilon);
            possibleMixPieceSegments.put(i, segmentsFromStartIdx);
        }
        double[][] best = new double[points.size()][];
        double angle = possibleMixPieceSegments.get(points.size() - 1).get(0).getAMax() - possibleMixPieceSegments.get(points.size() - 1).get(0).getAMax();
        best[points.size() - 1] = new double[]{1, Math.pow(angle, pow), 1};
        for (int i=points.size()-2; i>=0; i--) {
            Encoding.findBestWithAngle(i, possibleMixPieceSegments, best, pow);
        }

        int start = 0;
        int count = 0;
        while (start < points.size()) {
            segments.add(possibleMixPieceSegments.get(start).get((int) (best[start][2]-1)));
            start += (int) (best[start][2]) + 1;
            count++;
        }
        return segments;
    }


    private static void mergePerB(ArrayList<MixPieceSegment> segments, ArrayList<MixPieceSegment> mergedSegments, ArrayList<MixPieceSegment> unmergedSegments) {
        double aMinTemp = -Double.MAX_VALUE;
        double aMaxTemp = Double.MAX_VALUE;
        double b = Double.NaN;
        ArrayList<Long> timestamps = new ArrayList<>();

        segments.sort(Comparator.comparingDouble(MixPieceSegment::getB).thenComparingDouble(MixPieceSegment::getA));
        for (int i = 0; i < segments.size(); i++) {
            if (b != segments.get(i).getB()) {
                if (timestamps.size() == 1)
                    unmergedSegments.add(new MixPieceSegment(timestamps.get(0), aMinTemp, aMaxTemp, b));
                else {
                    for (Long timestamp : timestamps)
                        mergedSegments.add(new MixPieceSegment(timestamp, aMinTemp, aMaxTemp, b));
                }
                timestamps.clear();
                timestamps.add(segments.get(i).getInitTimestamp());
                aMinTemp = segments.get(i).getAMin();
                aMaxTemp = segments.get(i).getAMax();
                b = segments.get(i).getB();
                continue;
            }
            if (segments.get(i).getAMin() <= aMaxTemp && segments.get(i).getAMax() >= aMinTemp) {
                timestamps.add(segments.get(i).getInitTimestamp());
                aMinTemp = Math.max(aMinTemp, segments.get(i).getAMin());
                aMaxTemp = Math.min(aMaxTemp, segments.get(i).getAMax());
            } else {
                if (timestamps.size() == 1) unmergedSegments.add(segments.get(i - 1));
                else {
                    for (long timestamp : timestamps)
                        mergedSegments.add(new MixPieceSegment(timestamp, aMinTemp, aMaxTemp, b));
                }
                timestamps.clear();
                timestamps.add(segments.get(i).getInitTimestamp());
                aMinTemp = segments.get(i).getAMin();
                aMaxTemp = segments.get(i).getAMax();
            }
        }
        if (!timestamps.isEmpty()) {
            if (timestamps.size() == 1)
                unmergedSegments.add(new MixPieceSegment(timestamps.get(0), aMinTemp, aMaxTemp, b));
            else {
                for (long timestamp : timestamps)
                    mergedSegments.add(new MixPieceSegment(timestamp, aMinTemp, aMaxTemp, b));
            }
        }
    }

    private static void mergeAll(ArrayList<MixPieceSegment> segments, ArrayList<MixPieceSegment> mergedSegments, ArrayList<MixPieceSegment> unmergedSegments) {
        double aMinTemp = -Double.MAX_VALUE;
        double aMaxTemp = Double.MAX_VALUE;
        ArrayList<Double> bValues = new ArrayList<>();
        ArrayList<Long> timestamps = new ArrayList<>();

        segments.sort(Comparator.comparingDouble(MixPieceSegment::getAMin));
        for (int i = 0; i < segments.size(); i++) {
            if (segments.get(i).getAMin() <= aMaxTemp && segments.get(i).getAMax() >= aMinTemp) {
                timestamps.add(segments.get(i).getInitTimestamp());
                aMinTemp = Math.max(aMinTemp, segments.get(i).getAMin());
                aMaxTemp = Math.min(aMaxTemp, segments.get(i).getAMax());
                bValues.add(segments.get(i).getB());
            } else {
                if (timestamps.size() == 1) unmergedSegments.add(segments.get(i - 1));
                else {
                    for (int j = 0; j < timestamps.size(); j++)
                        mergedSegments.add(new MixPieceSegment(timestamps.get(j), aMinTemp, aMaxTemp, bValues.get(j)));
                }
                timestamps.clear();
                timestamps.add(segments.get(i).getInitTimestamp());
                aMinTemp = segments.get(i).getAMin();
                aMaxTemp = segments.get(i).getAMax();
                bValues.clear();
                bValues.add(segments.get(i).getB());
            }
        }
        if (!timestamps.isEmpty()) {
            if (timestamps.size() == 1)
                unmergedSegments.add(new MixPieceSegment(timestamps.get(0), aMinTemp, aMaxTemp, bValues.get(0)));
            else {
                for (int i = 0; i < timestamps.size(); i++)
                    mergedSegments.add(new MixPieceSegment(timestamps.get(i), aMinTemp, aMaxTemp, bValues.get(i)));
            }
        }
    }

    private static void merge(ArrayList<MixPieceSegment> segments) {
        perBSegments = new ArrayList<>();
        perASegments = new ArrayList<>();
        restSegments = new ArrayList<>();
        ArrayList<MixPieceSegment> temp = new ArrayList<>();

        mergePerB(segments, perBSegments, temp);
        if (!temp.isEmpty()) {
            mergeAll(temp, perASegments, restSegments);
        }
    }

    private static List<Point> toPoints() {
        ArrayList<MixPieceSegment> segments = new ArrayList<>();
        segments.addAll(perBSegments);
        segments.addAll(perASegments);
        segments.addAll(restSegments);
        segments.sort(Comparator.comparingLong(MixPieceSegment::getInitTimestamp));
        List<Point> points = new ArrayList<>();
        long currentTimeStamp = segments.get(0).getInitTimestamp();


        for (int i = 0; i < segments.size() - 1; i++) {
            while (currentTimeStamp < segments.get(i + 1).getInitTimestamp()) {
                points.add(new Point(currentTimeStamp, segments.get(i).getA() * (currentTimeStamp - segments.get(i).getInitTimestamp()) + segments.get(i).getB()));
                currentTimeStamp++;
            }
        }

        while (currentTimeStamp <= lastTimeStamp) {
            points.add(new Point(currentTimeStamp, segments.get(segments.size() - 1).getA() * (currentTimeStamp - segments.get(segments.size() - 1).getInitTimestamp()) + segments.get(segments.size() - 1).getB()));
            currentTimeStamp++;
        }

        return points;
    }


    private static void toByteArrayPerBSegmentsImproved(ArrayList<MixPieceSegment> segments, ByteArrayOutputStream outputStream, boolean variableByte) throws IOException {
        ByteArrayOutputStream tempStream = new ByteArrayOutputStream();
        List<Integer> groupSizes = new LinkedList<>();
        List<Integer> segmentSizes = new LinkedList<>();
        TreeMap<Integer, HashMap<Double, ArrayList<Long>>> input = new TreeMap<>();
        for (MixPieceSegment segment : segments) {
            double a = segment.getA();
            int b = (int) Math.round(segment.getB() / epsilon);
            long t = segment.getInitTimestamp();
            if (!input.containsKey(b)) input.put(b, new HashMap<>());
            if (!input.get(b).containsKey(a)) input.get(b).put(a, new ArrayList<>());
            input.get(b).get(a).add(t);
        }

        VariableByteEncoder.write(input.size(), tempStream);
//        if (variableByte) {
//            VariableByteEncoder.write(input.size(), outStream);
//        } else {
//            UIntEncoder.write(input.size(), outStream);
//        }
        if (input.isEmpty())
            return;
        int previousB = input.firstKey() - globalMinB;
        VariableByteEncoder.write(previousB, tempStream);
//        if (variableByte) {
//            VariableByteEncoder.write(previousB, outStream);
//        } else {
//            UIntEncoder.write(previousB, outStream);
//        }
        for (Map.Entry<Integer, HashMap<Double, ArrayList<Long>>> bSegments : input.entrySet()) {
            VariableByteEncoder.write(bSegments.getKey() - globalMinB - previousB, tempStream);
//            if (variableByte) {
//                VariableByteEncoder.write(bSegments.getKey() - globalMinB - previousB, outStream);
//            } else {
//                UIntEncoder.write(bSegments.getKey() - globalMinB - previousB, outStream);
//            }
            previousB = bSegments.getKey() - globalMinB;
            groupSizes.add(bSegments.getValue().size());
//            VariableByteEncoder.write(bSegments.getValue().size(), tempStream);
//            if (variableByte) {
//                VariableByteEncoder.write(bSegments.getValue().size(), tempStream);
//            } else {
//                UIntEncoder.write(bSegments.getValue().size(), tempStream);
//            }
            for (Map.Entry<Double, ArrayList<Long>> aSegment : bSegments.getValue().entrySet()) {
                FloatEncoder.write(aSegment.getKey().floatValue(), tempStream);
                Collections.sort(aSegment.getValue());
//                VariableByteEncoder.write(aSegment.getValue().size(), tempStream);
                segmentSizes.add(aSegment.getValue().size());
//                if (variableByte) {
//                    VariableByteEncoder.write(aSegment.getValue().size(), outStream);
//                } else {
//                    UIntEncoder.write(aSegment.getValue().size(), outStream);
//                }
                long previousTS = 0;
                for (Long timestamp : aSegment.getValue()) {
                    if (variableByte) {
                        VariableByteEncoder.write((int) (timestamp - previousTS), tempStream);
                    } else {
                        UIntEncoder.write((int) (timestamp - previousTS), tempStream);
                    }
                    previousTS = timestamp;
                }
            }
        }
        byte[] segmentSizesByteArray = getBinaryArrayFromList(segmentSizes);
        VariableByteEncoder.write(segmentSizesByteArray.length, outputStream);
        outputStream.write(segmentSizesByteArray);
        byte[] groupSizesByteArray = getBinaryArrayFromList(groupSizes);
        VariableByteEncoder.write(groupSizesByteArray.length, outputStream);
        outputStream.write(groupSizesByteArray);
        tempStream.writeTo(outputStream);
    }

    private static byte[] getBinaryArrayFromList(List<Integer> list) {
        int[] intArray = list.stream().mapToInt(i -> i).toArray();
        int[] compressed = new int[intArray.length];
        IntWrapper outputOffset = new IntWrapper(0);
        CODEC.compress(intArray, new IntWrapper(0), intArray.length, compressed, outputOffset);
        ByteBuffer byteBuffer = ByteBuffer.allocate(outputOffset.intValue() * 4);
        for (int i=0; i<outputOffset.intValue(); i++) {
            byteBuffer.putInt(compressed[i]); // Convert each int to 4 bytes
        }
        return byteBuffer.array();
    }

    private static void toByteArrayPerBSegments(ArrayList<MixPieceSegment> segments, ByteArrayOutputStream outStream, boolean variableByte) throws IOException {
        List<Integer> sizes = new LinkedList<>();
        TreeMap<Integer, HashMap<Double, ArrayList<Long>>> input = new TreeMap<>();
        for (MixPieceSegment segment : segments) {
            double a = segment.getA();
            int b = (int) Math.round(segment.getB() / epsilon);
            long t = segment.getInitTimestamp();
            if (!input.containsKey(b)) input.put(b, new HashMap<>());
            if (!input.get(b).containsKey(a)) input.get(b).put(a, new ArrayList<>());
            input.get(b).get(a).add(t);
        }

        VariableByteEncoder.write(input.size(), outStream);
//        if (variableByte) {
//            VariableByteEncoder.write(input.size(), outStream);
//        } else {
//            UIntEncoder.write(input.size(), outStream);
//        }
        if (input.isEmpty())
            return;
        int previousB = input.firstKey() - globalMinB;
        VariableByteEncoder.write(previousB, outStream);
//        if (variableByte) {
//            VariableByteEncoder.write(previousB, outStream);
//        } else {
//            UIntEncoder.write(previousB, outStream);
//        }
        for (Map.Entry<Integer, HashMap<Double, ArrayList<Long>>> bSegments : input.entrySet()) {
            VariableByteEncoder.write(bSegments.getKey() - globalMinB - previousB, outStream);
//            if (variableByte) {
//                VariableByteEncoder.write(bSegments.getKey() - globalMinB - previousB, outStream);
//            } else {
//                UIntEncoder.write(bSegments.getKey() - globalMinB - previousB, outStream);
//            }
            previousB = bSegments.getKey() - globalMinB;
            VariableByteEncoder.write(bSegments.getValue().size(), outStream);
//            if (variableByte) {
//                VariableByteEncoder.write(bSegments.getValue().size(), outStream);
//            } else {
//                UIntEncoder.write(bSegments.getValue().size(), outStream);
//            }
            for (Map.Entry<Double, ArrayList<Long>> aSegment : bSegments.getValue().entrySet()) {
                FloatEncoder.write(aSegment.getKey().floatValue(), outStream);
                Collections.sort(aSegment.getValue());
                VariableByteEncoder.write(aSegment.getValue().size(), outStream);
                sizes.add(aSegment.getValue().size());
//                if (variableByte) {
//                    VariableByteEncoder.write(aSegment.getValue().size(), outStream);
//                } else {
//                    UIntEncoder.write(aSegment.getValue().size(), outStream);
//                }
                long previousTS = 0;
                for (Long timestamp : aSegment.getValue()) {
                    if (variableByte) {
                        VariableByteEncoder.write((int) (timestamp - previousTS), outStream);
                    } else {
                        UIntEncoder.write((int) (timestamp - previousTS), outStream);
                    }
                    previousTS = timestamp;
                }
            }
        }
    }

    private static void toByteArrayPerASegmentsImproved(ArrayList<MixPieceSegment> segments, ByteArrayOutputStream outputStream, boolean variableByte) throws IOException {
        ByteArrayOutputStream tempStream = new ByteArrayOutputStream();
        List<Integer> betas = new LinkedList<>();
        TreeMap<Double, ArrayList<MixPieceSegment>> input = new TreeMap<>();
        for (MixPieceSegment segment : segments) {
            if (!input.containsKey(segment.getA())) input.put(segment.getA(), new ArrayList<>());
            input.get(segment.getA()).add(segment);
        }

        VariableByteEncoder.write(input.size(), tempStream);
//        if (variableByte) {
//            VariableByteEncoder.write(input.size(), outStream);
//        } else {
//            UIntEncoder.write(input.size(), outStream);
//        }
        for (Map.Entry<Double, ArrayList<MixPieceSegment>> aSegments : input.entrySet()) {
            FloatEncoder.write(aSegments.getKey().floatValue(), tempStream);
            VariableByteEncoder.write(aSegments.getValue().size(), tempStream);
//            if (variableByte) {
//                VariableByteEncoder.write(aSegments.getValue().size(), outStream);
//            } else {
//                UIntEncoder.write(aSegments.getValue().size(), outStream);
//            }
            aSegments.getValue().sort(Comparator.comparingDouble(MixPieceSegment::getB));
            int previousB = (int) Math.round(aSegments.getValue().get(0).getB() / epsilon) - globalMinB;
            VariableByteEncoder.write(previousB, tempStream);
//            if (variableByte) {
//                VariableByteEncoder.write(previousB, outStream);
//            } else {
//                UIntEncoder.write(previousB, outStream);
//            }
            for (MixPieceSegment segment : aSegments.getValue()) {
                int newB = (int) (Math.round(segment.getB() / epsilon) - globalMinB);
                betas.add(newB - previousB);
//                VariableByteEncoder.write(newB - previousB, tempStream);
//                if (variableByte) {
//                    VariableByteEncoder.write((int) (Math.round(segment.getB() / epsilon) - globalMinB - previousB), outStream);
//                } else {
//                    UIntEncoder.write((int) (Math.round(segment.getB() / epsilon) - globalMinB - previousB), outStream);
//                }
                previousB = newB;
                UIntEncoder.write(segment.getInitTimestamp(), tempStream);
            }
        }
        byte[] betasByteArray = getBinaryArrayFromList(betas);
        VariableByteEncoder.write(betasByteArray.length, outputStream);
        outputStream.write(betasByteArray);
        tempStream.writeTo(outputStream);
    }
    private static void toByteArrayPerASegments(ArrayList<MixPieceSegment> segments, ByteArrayOutputStream outStream, boolean variableByte) throws IOException {
        TreeMap<Double, ArrayList<MixPieceSegment>> input = new TreeMap<>();
        for (MixPieceSegment segment : segments) {
            if (!input.containsKey(segment.getA())) input.put(segment.getA(), new ArrayList<>());
            input.get(segment.getA()).add(segment);
        }

        VariableByteEncoder.write(input.size(), outStream);
//        if (variableByte) {
//            VariableByteEncoder.write(input.size(), outStream);
//        } else {
//            UIntEncoder.write(input.size(), outStream);
//        }
        for (Map.Entry<Double, ArrayList<MixPieceSegment>> aSegments : input.entrySet()) {
            FloatEncoder.write(aSegments.getKey().floatValue(), outStream);
            VariableByteEncoder.write(aSegments.getValue().size(), outStream);
//            if (variableByte) {
//                VariableByteEncoder.write(aSegments.getValue().size(), outStream);
//            } else {
//                UIntEncoder.write(aSegments.getValue().size(), outStream);
//            }
            aSegments.getValue().sort(Comparator.comparingDouble(MixPieceSegment::getB));
            int previousB = (int) Math.round(aSegments.getValue().get(0).getB() / epsilon) - globalMinB;
            VariableByteEncoder.write(previousB, outStream);
//            if (variableByte) {
//                VariableByteEncoder.write(previousB, outStream);
//            } else {
//                UIntEncoder.write(previousB, outStream);
//            }
            for (MixPieceSegment segment : aSegments.getValue()) {
                int newB = (int) (Math.round(segment.getB() / epsilon) - globalMinB);
                VariableByteEncoder.write(newB - previousB, outStream);
//                if (variableByte) {
//                    VariableByteEncoder.write((int) (Math.round(segment.getB() / epsilon) - globalMinB - previousB), outStream);
//                } else {
//                    UIntEncoder.write((int) (Math.round(segment.getB() / epsilon) - globalMinB - previousB), outStream);
//                }
                previousB = newB;
                UIntEncoder.write(segment.getInitTimestamp(), outStream);
            }
        }
    }

    private static void toByteArrayRestSegments(ArrayList<MixPieceSegment> segments, ByteArrayOutputStream outStream, boolean variableByte) throws IOException {
        List<Integer> betas = new LinkedList<>();
        VariableByteEncoder.write(segments.size(), outStream);
//        if (variableByte) {
//            VariableByteEncoder.write(segments.size(), outStream);
//        } else {
//            UIntEncoder.write(segments.size(), outStream);
//        }
        if (segments.isEmpty())
            return;
        segments.sort(Comparator.comparingDouble(MixPieceSegment::getB));
        int previousB = (int) Math.round(segments.get(0).getB() / epsilon) - globalMinB;
        VariableByteEncoder.write(previousB, outStream);
//        if (variableByte) {
//            VariableByteEncoder.write(previousB, outStream);
//        } else {
//            UIntEncoder.write(previousB, outStream);
//        }
        for (MixPieceSegment segment : segments) {
            betas.add((int) (Math.round(segment.getB() / epsilon) - globalMinB - previousB));
            VariableByteEncoder.write((int) (Math.round(segment.getB() / epsilon) - globalMinB - previousB), outStream);
//            if (variableByte) {
//                VariableByteEncoder.write((int) (Math.round(segment.getB() / epsilon) - globalMinB - previousB), outStream);
//            } else {
//                UIntEncoder.write((int) (Math.round(segment.getB() / epsilon) - globalMinB - previousB), outStream);
//            }
            previousB = (int) Math.round(segment.getB() / epsilon) - globalMinB;
            FloatEncoder.write((float) segment.getA(), outStream);
            UIntEncoder.write(segment.getInitTimestamp(), outStream);
        }
    }

    private static byte[] toByteArray(boolean variableByte, boolean zstd) {
        ByteArrayOutputStream outStream = new ByteArrayOutputStream();
        byte[] bytes = null;

        try {
            FloatEncoder.write((float) epsilon, outStream);
            VariableByteEncoder.write(globalMinB, outStream);
            toByteArrayPerBSegments(perBSegments, outStream, variableByte);
            toByteArrayPerASegments(perASegments, outStream, variableByte);
            toByteArrayRestSegments(restSegments, outStream, variableByte);

            if (variableByte) {
                VariableByteEncoder.write((int) lastTimeStamp, outStream);
            } else {
                UIntEncoder.write(lastTimeStamp, outStream);
            }

            if (zstd) {
                bytes = Zstd.compress(outStream.toByteArray());
            } else {
                bytes = outStream.toByteArray();
            }

            outStream.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        return bytes;
    }

    private static byte[] toByteArrayImproved(boolean variableByte, boolean zstd) {
        ByteArrayOutputStream outStream = new ByteArrayOutputStream();
        byte[] bytes = null;

        try {
            FloatEncoder.write((float) epsilon, outStream);
            VariableByteEncoder.write(globalMinB, outStream);
            toByteArrayPerBSegmentsImproved(perBSegments, outStream, variableByte);
            toByteArrayPerASegmentsImproved(perASegments, outStream, variableByte);
            toByteArrayRestSegments(restSegments, outStream, variableByte);

            if (variableByte) {
                VariableByteEncoder.write((int) lastTimeStamp, outStream);
            } else {
                UIntEncoder.write(lastTimeStamp, outStream);
            }

            if (zstd) {
                bytes = Zstd.compress(outStream.toByteArray());
            } else {
                bytes = outStream.toByteArray();
            }

            outStream.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        return bytes;
    }

    private static ArrayList<MixPieceSegment> readMergedPerBSegmentsImproved(ByteArrayInputStream inStream, boolean variableByte) throws IOException {
        ArrayList<MixPieceSegment> segments = new ArrayList<>();

        int[] recoveredSegmentSizes = recoverIntArray(inStream);
        int recoveredSegmentSizesIndex = 0;
        int[] recoveredGroupSizes = recoverIntArray(inStream);
        int recoveredGroupSizesIndex = 0;

        long numB = VariableByteEncoder.read(inStream);
//        long numB = variableByte ? VariableByteEncoder.read(inStream) : UIntEncoder.read(inStream);
        if (numB == 0)
            return segments;
        int previousB = VariableByteEncoder.read(inStream);
//        int previousB = variableByte ? VariableByteEncoder.read(inStream) : (int) UIntEncoder.read(inStream);
        for (int i = 0; i < numB; i++) {
            int b = VariableByteEncoder.read(inStream) + globalMinB + previousB;
//            int b = variableByte ? VariableByteEncoder.read(inStream) + globalMinB + previousB : (int) UIntEncoder.read(inStream) + globalMinB + previousB;
            previousB = b - globalMinB;
            int numA = recoveredGroupSizes[recoveredGroupSizesIndex++];
//            int numA = VariableByteEncoder.read(inStream);
//            System.out.println("Comparison: " + numA + " " + recoveredGroupSizes[recoveredGroupSizesIndex++]);
//            int numA = variableByte ? VariableByteEncoder.read(inStream) : (int) UIntEncoder.read(inStream);
            for (int j = 0; j < numA; j++) {
                float a = FloatEncoder.read(inStream);
                int numTimestamps = recoveredSegmentSizes[recoveredSegmentSizesIndex++];
//                int numTimestamps = VariableByteEncoder.read(inStream);
//                System.out.println("Comparison: " + numTimestamps + " " + recovered[noOfTimestampsIndex++]);
//                int numTimestamps = variableByte ? VariableByteEncoder.read(inStream) : (int) UIntEncoder.read(inStream);
                long timestamp = 0;
                for (int k = 0; k < numTimestamps; k++) {
                    timestamp += variableByte ? VariableByteEncoder.read(inStream) : UIntEncoder.read(inStream);
                    segments.add(new MixPieceSegment(timestamp, a, (float) (b * epsilon)));
                }
            }
        }
        return segments;
    }

    private static int[] recoverIntArray(ByteArrayInputStream inStream) throws IOException {
        int size = VariableByteEncoder.read(inStream);
        byte[] byteArray = new byte[size];
        inStream.read(byteArray, 0, size);
        int[] intArray = new int[byteArray.length / 4];
        ByteBuffer byteBuffer = ByteBuffer.wrap(byteArray);
        for (int i = 0; i < intArray.length; i++) {
            intArray[i] = byteBuffer.getInt();
        }
        int[] recovered = new int[size * 10];
        IntWrapper recoffset = new IntWrapper(0);
        CODEC.uncompress(intArray, new IntWrapper(0), intArray.length, recovered, recoffset);
        return recovered;
    }


    private static ArrayList<MixPieceSegment> readMergedPerBSegments(ByteArrayInputStream inStream, boolean variableByte) throws IOException {
        ArrayList<MixPieceSegment> segments = new ArrayList<>();
        long numB = VariableByteEncoder.read(inStream);
//        long numB = variableByte ? VariableByteEncoder.read(inStream) : UIntEncoder.read(inStream);
        if (numB == 0)
            return segments;
        int previousB = VariableByteEncoder.read(inStream);
//        int previousB = variableByte ? VariableByteEncoder.read(inStream) : (int) UIntEncoder.read(inStream);
        for (int i = 0; i < numB; i++) {
            int b = VariableByteEncoder.read(inStream) + globalMinB + previousB;
//            int b = variableByte ? VariableByteEncoder.read(inStream) + globalMinB + previousB : (int) UIntEncoder.read(inStream) + globalMinB + previousB;
            previousB = b - globalMinB;
            int numA = VariableByteEncoder.read(inStream);
//            int numA = variableByte ? VariableByteEncoder.read(inStream) : (int) UIntEncoder.read(inStream);
            for (int j = 0; j < numA; j++) {
                float a = FloatEncoder.read(inStream);
                int numTimestamps = VariableByteEncoder.read(inStream);
//                int numTimestamps = variableByte ? VariableByteEncoder.read(inStream) : (int) UIntEncoder.read(inStream);
                long timestamp = 0;
                for (int k = 0; k < numTimestamps; k++) {
                    timestamp += variableByte ? VariableByteEncoder.read(inStream) : UIntEncoder.read(inStream);
                    segments.add(new MixPieceSegment(timestamp, a, (float) (b * epsilon)));
                }
            }
        }
        return segments;
    }

    private static ArrayList<MixPieceSegment> readMergedPerASegmentsImproved(ByteArrayInputStream inStream, boolean variableByte) throws IOException {
        ArrayList<MixPieceSegment> segments = new ArrayList<>();
        int[] recoveredBetas = recoverIntArray(inStream);
        int recoveredBetasIndex = 0;

        int numA = VariableByteEncoder.read(inStream);
//        int numA = variableByte ? VariableByteEncoder.read(inStream) : (int) UIntEncoder.read(inStream);
        for (int i = 0; i < numA; i++) {
            float a = FloatEncoder.read(inStream);
            int numBT = VariableByteEncoder.read(inStream);
//            int numBT = variableByte ? VariableByteEncoder.read(inStream) : (int) UIntEncoder.read(inStream);
            int previousB = VariableByteEncoder.read(inStream);
//            int previousB = variableByte ? VariableByteEncoder.read(inStream) : (int) UIntEncoder.read(inStream);
            for (int j = 0; j < numBT; j++) {
                int b = recoveredBetas[recoveredBetasIndex++] + globalMinB + previousB;
//                int b = VariableByteEncoder.read(inStream) + globalMinB + previousB;
//                int b = variableByte ? VariableByteEncoder.read(inStream) + globalMinB + previousB : (int) (UIntEncoder.read(inStream) + globalMinB + previousB);
                previousB = b - globalMinB;
                long timestamp = UIntEncoder.read(inStream);
                segments.add(new MixPieceSegment(timestamp, a, (float) (b * epsilon)));
            }
        }

        return segments;
    }
    private static ArrayList<MixPieceSegment> readMergedPerASegments(ByteArrayInputStream inStream, boolean variableByte) throws IOException {
        ArrayList<MixPieceSegment> segments = new ArrayList<>();
        int numA = VariableByteEncoder.read(inStream);
//        int numA = variableByte ? VariableByteEncoder.read(inStream) : (int) UIntEncoder.read(inStream);
        for (int i = 0; i < numA; i++) {
            float a = FloatEncoder.read(inStream);
            int numBT = VariableByteEncoder.read(inStream);
//            int numBT = variableByte ? VariableByteEncoder.read(inStream) : (int) UIntEncoder.read(inStream);
            int previousB = VariableByteEncoder.read(inStream);
//            int previousB = variableByte ? VariableByteEncoder.read(inStream) : (int) UIntEncoder.read(inStream);
            for (int j = 0; j < numBT; j++) {
                int b = VariableByteEncoder.read(inStream) + globalMinB + previousB;
//                int b = variableByte ? VariableByteEncoder.read(inStream) + globalMinB + previousB : (int) (UIntEncoder.read(inStream) + globalMinB + previousB);
                previousB = b - globalMinB;
                long timestamp = UIntEncoder.read(inStream);
                segments.add(new MixPieceSegment(timestamp, a, (float) (b * epsilon)));
            }
        }

        return segments;
    }

    private static ArrayList<MixPieceSegment> readUnmerged(ByteArrayInputStream inStream, boolean variableByte) throws IOException {
        ArrayList<MixPieceSegment> segments = new ArrayList<>();
        int num = VariableByteEncoder.read(inStream);
//        int num = variableByte ? VariableByteEncoder.read(inStream) : (int) UIntEncoder.read(inStream);
        if (num == 0)
            return segments;
        int previousB = VariableByteEncoder.read(inStream);
//        int previousB = variableByte ? VariableByteEncoder.read(inStream) : (int) UIntEncoder.read(inStream);
        for (int i = 0; i < num; i++) {
            int b = VariableByteEncoder.read(inStream) + globalMinB + previousB;
//            int b = (variableByte ? VariableByteEncoder.read(inStream) : (int) UIntEncoder.read(inStream)) + globalMinB + previousB;
            previousB = b - globalMinB;
            float a = FloatEncoder.read(inStream);
            long timestamp = UIntEncoder.read(inStream);
            segments.add(new MixPieceSegment(timestamp, a, (float) (b * epsilon)));
        }

        return segments;
    }

    private static void readByteArray(byte[] input, boolean variableByte, boolean zstd) {
        byte[] binary;
        if (zstd) {
            binary = Zstd.decompress(input, input.length * 2); //TODO: How to know apriori original size?
        } else {
            binary = input;
        }
        ByteArrayInputStream inStream = new ByteArrayInputStream(binary);

        try {
            epsilon = FloatEncoder.read(inStream);
            globalMinB = VariableByteEncoder.read(inStream);

            perBSegments = readMergedPerBSegments(inStream, variableByte);
            perASegments = readMergedPerASegments(inStream, variableByte);
            restSegments = readUnmerged(inStream, variableByte);
            if (variableByte) {
                lastTimeStamp = VariableByteEncoder.read(inStream);
            } else {
                lastTimeStamp = UIntEncoder.read(inStream);
            }

            inStream.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static void readByteArrayImproved(byte[] input, boolean variableByte, boolean zstd) {
        byte[] binary;
        if (zstd) {
            binary = Zstd.decompress(input, input.length * 2); //TODO: How to know apriori original size?
        } else {
            binary = input;
        }
        ByteArrayInputStream inStream = new ByteArrayInputStream(binary);

        try {
            epsilon = FloatEncoder.read(inStream);
            globalMinB = VariableByteEncoder.read(inStream);

            perBSegments = readMergedPerBSegmentsImproved(inStream, variableByte);
            perASegments = readMergedPerASegmentsImproved(inStream, variableByte);
            restSegments = readUnmerged(inStream, variableByte);
            if (variableByte) {
                lastTimeStamp = VariableByteEncoder.read(inStream);
            } else {
                lastTimeStamp = UIntEncoder.read(inStream);
            }

            inStream.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}