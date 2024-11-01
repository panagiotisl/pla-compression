package org.pla.compression.encodings;

import com.github.benmanes.caffeine.cache.Caffeine;
import com.github.benmanes.caffeine.cache.LoadingCache;
import com.github.luben.zstd.Zstd;
import org.pla.compression.util.Encoding.FloatEncoder;
import org.pla.compression.util.Encoding.UIntEncoder;
import org.pla.compression.util.Encoding.VariableByteEncoder;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.*;

public class Encoding {
    private static final double POWER = 0.1;
    private List<MixPieceSegment> segments;
    private byte[] byteArray;

    private List<MixPieceSegment> bestMixPieceSegments;

    private int bestSize = Integer.MAX_VALUE;
    private double epsilon;
    private long lastTimeStamp;

    private int previousNoOfGroups = 0;
    private List<Point> points;

    public Encoding(List<Point> points, double epsilon, int mode, double pow) throws IOException {
        if (points.isEmpty()) throw new IOException();

        this.points = points;
        this.epsilon = epsilon;
        this.lastTimeStamp = points.get(points.size() - 1).getTimestamp();
        this.segments = compress(points, mode, pow);
    }

    public Encoding(byte[] bytes, boolean variableByte, boolean zstd) throws IOException {
        readByteArray(bytes, variableByte, zstd);
    }

    private static double quantization(double value, double epsilon) {
        return Math.round(value / epsilon) * epsilon;
    }


    public static double ceil(double value, double epsilon) {
        return Math.ceil(value / epsilon) * epsilon;
    }
    public static double floor(double value, double epsilon) {
        return Math.floor(value / epsilon) * epsilon;
    }

    public List<MixPieceSegment> getSegments() {
        return segments;
    }

    public static List<MixPieceSegment> createMixPieceSegmentsFromStartIdx(int startIdx, List<Point> points, double epsilon) {
        List<MixPieceSegment> segments = new LinkedList<>();

        long initTimestamp = points.get(startIdx).getTimestamp();
        double b1 = floor(points.get(startIdx).getValue(), epsilon);
        double b2 = ceil(points.get(startIdx).getValue(), epsilon);
        double b = quantization(points.get(startIdx).getValue(), epsilon);
        boolean floor = true;
        boolean ceil = true;
        int length = 0;
        if (startIdx + 1 == points.size()) {
            segments.add(new MixPieceSegment(initTimestamp, -Double.MAX_VALUE, Double.MAX_VALUE, b));
            return segments;
        }
        double aMax1 = ((points.get(startIdx + 1).getValue() + epsilon) - b1) / (points.get(startIdx + 1).getTimestamp() - initTimestamp);
        double aMin1 = ((points.get(startIdx + 1).getValue() - epsilon) - b1) / (points.get(startIdx + 1).getTimestamp() - initTimestamp);
        double aMax2 = ((points.get(startIdx + 1).getValue() + epsilon) - b2) / (points.get(startIdx + 1).getTimestamp() - initTimestamp);
        double aMin2 = ((points.get(startIdx + 1).getValue() - epsilon) - b2) / (points.get(startIdx + 1).getTimestamp() - initTimestamp);

        if (startIdx + 2 == points.size()) {
            segments.add(new MixPieceSegment(initTimestamp, aMin1, aMax1, b1));
            return segments;
        }

        for (int idx = startIdx + 2; idx < points.size(); idx++) {
            double upValue = points.get(idx).getValue() + epsilon;
            double downValue = points.get(idx).getValue() - epsilon;

            double upLim1 = aMax1 * (points.get(idx).getTimestamp() - initTimestamp) + b1;
            double downLim1 = aMin1 * (points.get(idx).getTimestamp() - initTimestamp) + b1;
            double upLim2 = aMax2 * (points.get(idx).getTimestamp() - initTimestamp) + b2;
            double downLim2 = aMin2 * (points.get(idx).getTimestamp() - initTimestamp) + b2;
            if ((downValue > upLim1 || upValue < downLim1)) {
                floor = false;
            }
            if ((downValue > upLim2 || upValue < downLim2)) {
                ceil = false;
            }
            if (floor) length++;
            if (ceil) length--;
            if (length > 0) {
                segments.add(new MixPieceSegment(initTimestamp, aMin1, aMax1, b1));
            } else if (length < 0){
                segments.add(new MixPieceSegment(initTimestamp, aMin2, aMax2, b2));
            } else if (aMax1 - aMin1 > aMax2 - aMin2) {
                segments.add(new MixPieceSegment(initTimestamp, aMin1, aMax1, b1));
            } else {
                segments.add(new MixPieceSegment(initTimestamp, aMin2, aMax2, b2));
            }
            if (!floor && !ceil) {
                return segments;
            }

            if (upValue < upLim1)
                aMax1 = Math.max((upValue - b1) / (points.get(idx).getTimestamp() - initTimestamp), aMin1);
            if (downValue > downLim1)
                aMin1 = Math.min((downValue - b1) / (points.get(idx).getTimestamp() - initTimestamp), aMax1);
            if (upValue < upLim2)
                aMax2 = Math.max((upValue - b2) / (points.get(idx).getTimestamp() - initTimestamp), aMin2);
            if (downValue > downLim2)
                aMin2 = Math.min((downValue - b2) / (points.get(idx).getTimestamp() - initTimestamp), aMax2);
        }
        if (length > 0) {
            segments.add(new MixPieceSegment(initTimestamp, aMin1, aMax1, b1));
        } else if (length < 0){
            segments.add(new MixPieceSegment(initTimestamp, aMin2, aMax2, b2));
        } else if (aMax1 - aMin1 > aMax2 - aMin2) {
            segments.add(new MixPieceSegment(initTimestamp, aMin1, aMax1, b1));
        } else {
            segments.add(new MixPieceSegment(initTimestamp, aMin2, aMax2, b2));
        }

        return segments;
    }

    public static int findReachFrom(int startIdx, List<Point> points, double epsilon) {
        if (startIdx + 1 == points.size()) {
            return 1;
        }
        if (startIdx + 2 == points.size()) {
            return 1;
        }
        int segments = 0;
        boolean floor = true;
        boolean ceil = true;
        long initTimestamp = points.get(startIdx).getTimestamp();
        double b1 = floor(points.get(startIdx).getValue(), epsilon);
        double b2 = ceil(points.get(startIdx).getValue(), epsilon);
        double aMax1 = ((points.get(startIdx + 1).getValue() + epsilon) - b1) / (points.get(startIdx + 1).getTimestamp() - initTimestamp);
        double aMin1 = ((points.get(startIdx + 1).getValue() - epsilon) - b1) / (points.get(startIdx + 1).getTimestamp() - initTimestamp);
        double aMax2 = ((points.get(startIdx + 1).getValue() + epsilon) - b2) / (points.get(startIdx + 1).getTimestamp() - initTimestamp);
        double aMin2 = ((points.get(startIdx + 1).getValue() - epsilon) - b2) / (points.get(startIdx + 1).getTimestamp() - initTimestamp);

        for (int idx = startIdx + 2; idx < points.size(); idx++) {
            double upValue = points.get(idx).getValue() + epsilon;
            double downValue = points.get(idx).getValue() - epsilon;

            double upLim1 = aMax1 * (points.get(idx).getTimestamp() - initTimestamp) + b1;
            double downLim1 = aMin1 * (points.get(idx).getTimestamp() - initTimestamp) + b1;
            double upLim2 = aMax2 * (points.get(idx).getTimestamp() - initTimestamp) + b2;
            double downLim2 = aMin2 * (points.get(idx).getTimestamp() - initTimestamp) + b2;
            segments++;
            if ((downValue > upLim1 || upValue < downLim1)) {
                floor = false;
            }
            if ((downValue > upLim2 || upValue < downLim2)) {
                ceil = false;
            }
            if (!floor && !ceil) {
                return segments;
            }
            if (upValue < upLim1)
                aMax1 = Math.max((upValue - b1) / (points.get(idx).getTimestamp() - initTimestamp), aMin1);
            if (downValue > downLim1)
                aMin1 = Math.min((downValue - b1) / (points.get(idx).getTimestamp() - initTimestamp), aMax1);
            if (upValue < upLim2)
                aMax2 = Math.max((upValue - b2) / (points.get(idx).getTimestamp() - initTimestamp), aMin2);
            if (downValue > downLim2)
                aMin2 = Math.min((downValue - b2) / (points.get(idx).getTimestamp() - initTimestamp), aMax2);
        }
        segments++;
        return segments;
    }


    private State binarySearch(List<State> states, int previousNoOfGroups, long initTimestamp, double b, List<MixPieceSegment> segments) {
        List<MixPieceSegment> oldMixPieceSegments = mergePerB(segments);
        TreeMap<Double, HashMap<Double, List<Long>>> tree = new TreeMap<>();
        for (MixPieceSegment segment : oldMixPieceSegments) {
            double a = segment.getA();
            double bi = segment.getB();
            long t = segment.getInitTimestamp();
            if (!tree.containsKey(bi)) tree.put(bi, new HashMap<>());
            if (!tree.get(bi).containsKey(a)) tree.get(bi).put(a, new ArrayList<>());
            tree.get(bi).get(a).add(t);
        }

        int high = states.size() - 1;
        int low = high / 2;
        State found = states.get(high);

        while (low <= high) {
            int mid = low  + ((high - low) / 2);
            State state = states.get(mid);
            if (findNumberOfGroups(initTimestamp, state.aMax, state.aMin, b, segments, previousNoOfGroups, oldMixPieceSegments, tree)) {
                break;
            } else {
                found = state;
                low = mid + 1;
            }
        }
        return found;
    }

    private class State {

        int idx;

        double aMin;

        double aMax;

        public State(final int idx, final double aMin, final double aMax) {
            this.idx = idx;
            this.aMin = aMin;
            this.aMax = aMax;
        }

        @Override
        public String toString() {
            return String.format("%d, %f, %f, %f", idx, aMin, aMax, aMax - aMin);
        }
    }
    private int findNumberOfGroups(long initTimestamp, double aMax, double aMin, double b, List<MixPieceSegment> segments) {
        List<MixPieceSegment> tempMixPieceSegments = new ArrayList<>(segments);
        tempMixPieceSegments.add(new MixPieceSegment(initTimestamp, aMin, aMax, b));
        tempMixPieceSegments = mergePerB(tempMixPieceSegments);

        TreeMap<Double, HashMap<Double, ArrayList<Long>>> input = new TreeMap<>();
        for (MixPieceSegment segment : tempMixPieceSegments) {
            double a = segment.getA();
            double bi = segment.getB();
            long t = segment.getInitTimestamp();
            if (!input.containsKey(bi)) input.put(bi, new HashMap<>());
            if (!input.get(bi).containsKey(a)) input.get(bi).put(a, new ArrayList<>());
            input.get(bi).get(a).add(t);
        }
        int groups = 0;
        for (Map.Entry<Double, HashMap<Double, ArrayList<Long>>> bMixPieceSegments : input.entrySet())
            groups += bMixPieceSegments.getValue().size();
        return groups;
    }

    private boolean findNumberOfGroups(long initTimestamp, double aMax, double aMin, double b, List<MixPieceSegment> segments,
                                       int previousNoOfGroups, List<MixPieceSegment> oldMixPieceSegments,
                                       TreeMap<Double, HashMap<Double, List<Long>>> tree) {
        double alpha = oldMixPieceSegments.stream().filter(s -> s.getB() == b && s.getAMin() <= aMax && s.getAMax() >= aMin).mapToDouble(s -> s.getAMax() - s.getAMin()).max().orElse(-1.0);

        List<MixPieceSegment> tempMixPieceSegments = new ArrayList<>(segments);
        tempMixPieceSegments.add(new MixPieceSegment(initTimestamp, aMin, aMax, b));
        tempMixPieceSegments = mergePerB(tempMixPieceSegments);
        double value = 0;
        TreeMap<Double, HashMap<Double, List<Long>>> input = new TreeMap<>();
        for (MixPieceSegment segment : tempMixPieceSegments) {
            double a = segment.getA();
            double bi = segment.getB();
            long t = segment.getInitTimestamp();
            if (initTimestamp == t) {
//                System.out.println(String.format("A: %f, A_new: %f, B: %f, A/A_new: %f", alpha, (segment.getAMax() - segment.getAMin()), (aMax - aMin),
//                        alpha / (segment.getAMax() - segment.getAMin())));
                value = alpha / (segment.getAMax() - segment.getAMin());
            }
            if (!input.containsKey(bi)) input.put(bi, new HashMap<>());
            if (!input.get(bi).containsKey(a)) input.get(bi).put(a, new ArrayList<>());
            input.get(bi).get(a).add(t);
        }
        int groups = 0;
        for (Map.Entry<Double, HashMap<Double, List<Long>>> bMixPieceSegments : input.entrySet())
            groups += bMixPieceSegments.getValue().size();

        return groups > previousNoOfGroups;// || value > 10;
    }

    LoadingCache<Integer, List<MixPieceSegment>> cache = Caffeine.newBuilder()
            .maximumSize(5000)
            .build(this::createMixPieceSegmentsFromStartIdx);
//    Map<Integer, List<SimPieceMixPieceSegment>> cache = new HashMap<>();

    private List<MixPieceSegment> createMixPieceSegmentsFromStartIdx(int idx) {
        return createMixPieceSegmentsFromStartIdx(idx, this.points, epsilon);
    }

    private List<MixPieceSegment> compress(List<Point> points) {
        return compress(points, 0, 0);
    }

    private List<MixPieceSegment> compress(List<Point> points, int mode, double pow) {
        Map<Integer, List<MixPieceSegment>> possibleMixPieceSegments = new TreeMap<>();
        switch (mode) {
            case 0:
                for (int i = 0; i <= points.size() - 1; i++) {
                    List<MixPieceSegment> segmentsFromStartIdx = createMixPieceSegmentsFromStartIdx(i, points, epsilon);
//            System.out.println(String.format("StartIdx: %d: %d", i, segmentsFromStartIdx.size()));
                    possibleMixPieceSegments.put(i, segmentsFromStartIdx);
                }
                this.segments = new ArrayList<>();
                double[][] best = new double[points.size()][];
                double angle = possibleMixPieceSegments.get(points.size() - 1).get(0).getAMax() - possibleMixPieceSegments.get(points.size() - 1).get(0).getAMax();
                best[points.size() - 1] = new double[]{1, Math.pow(angle, pow), 1};
                for (int i=points.size()-2; i>=0; i--) {
                    findBestWithAngle(i, possibleMixPieceSegments, best, pow);
                }

                int start = 0;
                int count = 0;
                while (start < points.size()) {
                    this.segments.add(possibleMixPieceSegments.get(start).get((int) (best[start][2]-1)));
                    start += (int) (best[start][2]) + 1;
                    count++;
                }

//                System.out.println("Best: " + best[0][0] + " Used: " + this.segments.size());

                this.segments = mergePerB(this.segments);
                break;
            case 1:
                List<MixPieceSegment> segments = new ArrayList<>();
                int currentIdx = 0;
                while (currentIdx < points.size()) currentIdx = createMixPieceSegment(currentIdx, points, segments);
                this.segments = mergePerB(segments);
                break;
            case 2:
                this.segments = new ArrayList<>();
                int startIdx = 0;
                while (startIdx < points.size()) {
                    startIdx = addSegment(startIdx, pow, this.points, epsilon, this.segments);
                }
                this.segments = mergePerB(this.segments);
                break;
            case 3:
                for (int i = 0; i <= points.size() - 1; i++) {
                    List<MixPieceSegment> segmentsFromStartIdx = createMixPieceSegmentsFromStartIdx(i, points, epsilon);
                    possibleMixPieceSegments.put(i, segmentsFromStartIdx);
                }
                this.segments = new ArrayList<>();
                int startIdxG = 0;
                while (startIdxG < points.size()) {
                    startIdxG = addSegment(startIdxG, possibleMixPieceSegments, pow);
                }
                this.segments = mergePerB(this.segments);
                break;
            case 4:
                for (int i = 0; i <= points.size() - 1; i++) {
                    List<MixPieceSegment> segmentsFromStartIdx = createMixPieceSegmentsFromStartIdx(i, points, epsilon);
                    possibleMixPieceSegments.put(i, segmentsFromStartIdx);
                }
                this.segments = new ArrayList<>();
                startIdx = 0;
                while (startIdx < points.size()) {
                    startIdx = addMixPieceSegment2(startIdx, pow);
                }
                this.segments = mergePerB(this.segments);
                break;
            case 5:
                List<MixPieceSegment> segmentsBoth = new ArrayList<>();
                int currentIdxBoth = 0;
                while (currentIdxBoth < points.size()) currentIdxBoth = createMixPieceSegmentBoth(currentIdxBoth, points, segmentsBoth);
                this.segments = mergePerB(segmentsBoth);
                break;
        }
        return segments;
    }


    private int createMixPieceSegment(int startIdx, List<Point> points, List<MixPieceSegment> segments) {
        long initTimestamp = points.get(startIdx).getTimestamp();
        double b = floor(points.get(startIdx).getValue(), epsilon);
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

    private int createMixPieceSegmentBoth(int startIdx, List<Point> points, List<MixPieceSegment> segments) {
        long initTimestamp = points.get(startIdx).getTimestamp();
        double b1 = floor(points.get(startIdx).getValue(), epsilon);
        double b2 = ceil(points.get(startIdx).getValue(), epsilon);
        int count = 0;
        boolean floor = true;
        boolean ceil = true;
        double b = quantization(points.get(startIdx).getValue(), epsilon);
        if (startIdx + 1 == points.size()) {
            segments.add(new MixPieceSegment(initTimestamp, -Double.MAX_VALUE, Double.MAX_VALUE, b));
            return startIdx + 1;
        }
        double aMax1 = ((points.get(startIdx + 1).getValue() + epsilon) - b1) / (points.get(startIdx + 1).getTimestamp() - initTimestamp);
        double aMin1 = ((points.get(startIdx + 1).getValue() - epsilon) - b1) / (points.get(startIdx + 1).getTimestamp() - initTimestamp);
        double aMax2 = ((points.get(startIdx + 1).getValue() + epsilon) - b2) / (points.get(startIdx + 1).getTimestamp() - initTimestamp);
        double aMin2 = ((points.get(startIdx + 1).getValue() - epsilon) - b2) / (points.get(startIdx + 1).getTimestamp() - initTimestamp);
        if (startIdx + 2 == points.size()) {
            segments.add(new MixPieceSegment(initTimestamp, aMin1, aMax1, b));
            return startIdx + 2;
        }

        for (int idx = startIdx + 2; idx < points.size(); idx++) {
            double upValue = points.get(idx).getValue() + epsilon;
            double downValue = points.get(idx).getValue() - epsilon;

            double upLim1 = aMax1 * (points.get(idx).getTimestamp() - initTimestamp) + b1;
            double downLim1 = aMin1 * (points.get(idx).getTimestamp() - initTimestamp) + b1;
            double upLim2 = aMax2 * (points.get(idx).getTimestamp() - initTimestamp) + b2;
            double downLim2 = aMin2 * (points.get(idx).getTimestamp() - initTimestamp) + b2;
            if ((downValue > upLim1 || upValue < downLim1)) {
                floor = false;
            }
            if ((downValue > upLim2 || upValue < downLim2)) {
                ceil = false;
            }
            if (floor) count++;
            if (ceil) count--;
            if (!floor && !ceil) {
                if (count > 0) {
                    segments.add(new MixPieceSegment(initTimestamp, aMin1, aMax1, b1));
                } else {
                    segments.add(new MixPieceSegment(initTimestamp, aMin2, aMax2, b2));
                }
                return idx;
            }
            if (upValue < upLim1)
                aMax1 = Math.max((upValue - b1) / (points.get(idx).getTimestamp() - initTimestamp), aMin1);
            if (downValue > downLim1)
                aMin1 = Math.min((downValue - b1) / (points.get(idx).getTimestamp() - initTimestamp), aMax1);
            if (upValue < upLim2)
                aMax2 = Math.max((upValue - b2) / (points.get(idx).getTimestamp() - initTimestamp), aMin2);
            if (downValue > downLim2)
                aMin2 = Math.min((downValue - b2) / (points.get(idx).getTimestamp() - initTimestamp), aMax2);
        }
        if (count > 0) {
            segments.add(new MixPieceSegment(initTimestamp, aMin1, aMax1, b1));
        } else {
            segments.add(new MixPieceSegment(initTimestamp, aMin2, aMax2, b2));
        }

        return points.size();
    }



    private int createMixPieceSegmentSimPath(int startIdx, List<Point> points, List<MixPieceSegment> segments) {
//        System.out.println("Starting segment " + segments.size());
        long initTimestamp = points.get(startIdx).getTimestamp();
        double b = quantization(points.get(startIdx).getValue(), epsilon);
        if (startIdx + 1 == points.size()) {
            segments.add(new MixPieceSegment(initTimestamp, -Double.MAX_VALUE, Double.MAX_VALUE, b));
            return startIdx + 1;
        }
        double aMax = ((points.get(startIdx + 1).getValue() + epsilon) - b) / (points.get(startIdx + 1).getTimestamp() - initTimestamp);
        double aMin = ((points.get(startIdx + 1).getValue() - epsilon) - b) / (points.get(startIdx + 1).getTimestamp() - initTimestamp);
//      System.out.println(aMax-aMin);
        double diff = aMax-aMin;
        if (startIdx + 2 == points.size()) {
            segments.add(new MixPieceSegment(initTimestamp, aMin, aMax, b));
            return startIdx + 2;
        }

        int groups, startGroups = 0;
        List<State> states = new ArrayList<>();
        for (int idx = startIdx + 2; idx < points.size(); idx++) {
//            if (idx - startIdx == 3) {
//                groups = findNumberOfGroups(initTimestamp, aMax, aMin, b, segments);
//                startGroups = groups;
//            }
            double upValue = points.get(idx).getValue() + epsilon;
            double downValue = points.get(idx).getValue() - epsilon;

            double upLim = aMax * (points.get(idx).getTimestamp() - initTimestamp) + b;
            double downLim = aMin * (points.get(idx).getTimestamp() - initTimestamp) + b;

            double aMaxTemp = aMax, aMinTemp = aMin;
            if (upValue < upLim)
                aMaxTemp = Math.max((upValue - b) / (points.get(idx).getTimestamp() - initTimestamp), aMin);
            if (downValue > downLim)
                aMinTemp = Math.min((downValue - b) / (points.get(idx).getTimestamp() - initTimestamp), aMax);
            states.add(new State(idx, aMin, aMax));
            if ((downValue > upLim || upValue < downLim)) {
                int newGroups = findNumberOfGroups(initTimestamp, aMax, aMin, b, segments);
                if (newGroups > previousNoOfGroups) {
//                    System.out.println(String.format("MixPieceSegments: %d, Groups: %d, Size: %d", segments.size(), newGroups, (idx - startIdx)));
                    State state = binarySearch(states, previousNoOfGroups, initTimestamp, b, segments);
//                    if (idx > state.idx) {
//                        System.out.println(String.format("%d - %d - %f", (idx - startIdx), (state.idx - startIdx) , ((double) state.idx - startIdx) / (idx - startIdx)));
//                    }
                    previousNoOfGroups = newGroups;
                    segments.add(new MixPieceSegment(initTimestamp, state.aMin, state.aMax, b));
                    return state.idx;
                }
                previousNoOfGroups = newGroups;
                segments.add(new MixPieceSegment(initTimestamp, aMin, aMax, b));
                return idx;
            }
            diff = aMax-aMin;
            if (upValue < upLim)
                aMax = Math.max((upValue - b) / (points.get(idx).getTimestamp() - initTimestamp), aMin);
            if (downValue > downLim)
                aMin = Math.min((downValue - b) / (points.get(idx).getTimestamp() - initTimestamp), aMax);
            diff = diff - (aMax-aMin);
        }
        segments.add(new MixPieceSegment(initTimestamp, aMin, aMax, b));

        return points.size();
    }

    private int findBest(int start, Map<Integer, List<MixPieceSegment>> possibleMixPieceSegments, int[][] best) {
//        System.out.println("CALLED: " + start + " - " + best.length);
        if (start >= possibleMixPieceSegments.size()) {
            return 0;
        }
        if (best[start] != null) {
            return best[start][0];
        }
        else {
            int bestResult = Integer.MAX_VALUE;
            int bestIndex = 0;
            for (int i=1; i<= possibleMixPieceSegments.get(start).size(); i++) {
                int result = 1 + findBest(start + i + 1, possibleMixPieceSegments, best);
                if (result < bestResult) {
                    bestResult = result;
                    bestIndex = i;
                }
            }
            best[start] = new int[]{bestResult, bestIndex};
//            System.out.println("Best for: " + start + " is: " + bestResult);
            return bestResult;
        }
    }



    public static double[] findBestWithAngle(int start, Map<Integer, List<MixPieceSegment>> possibleMixPieceSegments, double[][] best, double pow) {
        if (start >= possibleMixPieceSegments.size()) {
            return new double[] {0, 0};
        }
        if (best[start] != null) {
            return new double[] {best[start][0], best[start][1]};
        }
        else {
            double bestResult = Double.MIN_VALUE;
            double bestAngle = 0;
            double bestN = 0;
            int bestIndex = 0;
            for (int i=1; i<= possibleMixPieceSegments.get(start).size(); i++) {
                double[] result = findBestWithAngle(start + i + 1, possibleMixPieceSegments, best, pow);
                double n = result[0] + 1;
                MixPieceSegment segment = possibleMixPieceSegments.get(start).get(i-1);
                double angle = result[1] + Math.pow((segment.getAMax() - segment.getAMin()), pow);
                double cost = angle / (n * n);
//                double cost = 1.0 / (n * n);
//                double cost = 1.0 / n;
                if (cost > bestResult) {
                    bestResult = cost;
                    bestAngle = angle;
                    bestN = n;
                    bestIndex = i;
                }
            }
            best[start] = new double[]{bestN, bestAngle, bestIndex};
            return new double[] {bestN, bestAngle};
        }
    }

    public static int addSegment(int startIdx, double pow, List<Point> points, double epsilon, List<MixPieceSegment> segments) {
//        int firstMixPieceSegments = cache.get(startIdx).size();
        List<MixPieceSegment> firstMixPieceSegments = createMixPieceSegmentsFromStartIdx(startIdx, points, epsilon);
        int index = 0;
        if (startIdx + 2 < points.size()) {
//            int best = 2 + cache.get(startIdx + 2).size();
            int best = 2 + findReachFrom(startIdx + 2, points, epsilon);
            for (int i = (int) Math.pow(firstMixPieceSegments.size(), pow); i<firstMixPieceSegments.size() && startIdx + i + 2 < points.size(); i++) {
//                int reach = 2 + i + cache.get(startIdx + i + 2).size();
                int reach = 2 + i + findReachFrom(startIdx + i + 2, points, epsilon);
                if (reach > best) {
                    best = reach;
                    index = i;
                }
            }
        }
//        this.segments.add(cache.get(startIdx).get(index));
        segments.add(firstMixPieceSegments.get(index));
        return startIdx + index + 1 + 1;
    }


    private int addMixPieceSegment2(int startIdx, double pow) {


//        while (sequence not finished) {
//            i++;
//            start segment i;
//            extend segment i as far as allowed by error bound ε;
//            b_last = end boundary of segment i;
//            b_prev = end boundary of segment i-1;
//            b_start = start boundary of segment i-1;
//            sum_angles = (α_{i-1})^p + (α_i)^p;
//            best_prev = b_prev;
//            while (segment(b_prev--, b_last) allowed by ε) {
//                new_last_seg = segment(b_prev, b_last);
//                new_prev_seg = segment(b_start, b_prev);
//                if α(new_last_seg)^p + α(new_prev_seg)^p > sum_angles {
//                    sum_angles = α(new_last_seg)^p + α(new_prev_seg)^p;
//                    best_prev = b_prev;
//                }
//
//            }

        List<MixPieceSegment> firstMixPieceSegments = createMixPieceSegmentsFromStartIdx(startIdx);
        MixPieceSegment currMixPieceSegment = firstMixPieceSegments.get(firstMixPieceSegments.size() - 1);
        if (!this.segments.isEmpty()) {
            MixPieceSegment prevMixPieceSegment = this.segments.get(this.segments.size() - 1);
            List<MixPieceSegment> startPrevMixPieceSegments = createMixPieceSegmentsFromStartIdx((int) prevMixPieceSegment.getInitTimestamp());
            int bLast = startIdx + firstMixPieceSegments.size();
            int bPrev = startIdx - 1;
            int bStart = (int) prevMixPieceSegment.getInitTimestamp();
            double sumAngles = Math.pow(prevMixPieceSegment.getAMax() - prevMixPieceSegment.getAMin(), pow) + Math.pow(currMixPieceSegment.getAMax() - currMixPieceSegment.getAMin(), pow);
            int bestPrev = bPrev;
            while (bPrev >= 0) {
                List<MixPieceSegment> prevMixPieceSegments = createMixPieceSegmentsFromStartIdx(bPrev);
                if (bPrev + prevMixPieceSegments.size() + 1 < bLast) {
                    break;
                }
                if (bLast - bPrev - 1 >= prevMixPieceSegments.size() || bPrev - bStart - 1 - 1 < 1) {
                    bPrev--;
                    continue;
                }
                MixPieceSegment newLastMixPieceSegment = prevMixPieceSegments.get(bLast - bPrev - 1);
                MixPieceSegment newPrevMixPieceSegment = startPrevMixPieceSegments.get(bPrev - bStart - 1 - 1);
                double newSum = Math.pow(newLastMixPieceSegment.getAMax() - newLastMixPieceSegment.getAMin(), pow) + Math.pow(newPrevMixPieceSegment.getAMax() - newPrevMixPieceSegment.getAMin(), pow);
                if (newSum > sumAngles) {
                    sumAngles = newSum;
                    bestPrev = bPrev;
                    prevMixPieceSegment = newPrevMixPieceSegment;
                    currMixPieceSegment = newLastMixPieceSegment;
                }
                bPrev--;
            }
            this.segments.remove(this.segments.size()-1);
            this.segments.add(prevMixPieceSegment);
        }
        this.segments.add(currMixPieceSegment);
        return startIdx + firstMixPieceSegments.size() + 1;
    }


    private int addSegment(int startIdx, Map<Integer, List<MixPieceSegment>> possibleMixPieceSegments, double pow) {
//        System.out.println("Getting: " + startIdx);
//        System.out.println("Size: " + possibleMixPieceSegments.get(startIdx).size());
        int firstMixPieceSegments = possibleMixPieceSegments.get(startIdx).size();
        int index = 0;
        if (startIdx + 2 < possibleMixPieceSegments.size()) {
            double alpha1 = possibleMixPieceSegments.get(startIdx).get(0).getAMax() - possibleMixPieceSegments.get(startIdx).get(0).getAMax();
            double bestAlpha = Math.pow(alpha1, pow) + Math.pow(
                    possibleMixPieceSegments.get(startIdx + 2).get(possibleMixPieceSegments.get(startIdx + 2).size() -1).getAMax() -
                            possibleMixPieceSegments.get(startIdx + 2).get(possibleMixPieceSegments.get(startIdx + 2).size() -1).getAMin(), pow);
            int best = 2 + possibleMixPieceSegments.get(startIdx + 2).size();
//            System.out.println(0 + ": " + startIdx + " " + best);
            for (int i=1; i<firstMixPieceSegments && startIdx + i + 2 < possibleMixPieceSegments.size(); i++) {
                double alphaTotal = Math.pow(alpha1, pow) + Math.pow(
                        possibleMixPieceSegments.get(startIdx + i + 2).get(possibleMixPieceSegments.get(startIdx + i + 2).size() -1).getAMax() -
                                possibleMixPieceSegments.get(startIdx + i + 2).get(possibleMixPieceSegments.get(startIdx + i + 2).size() -1).getAMin(), pow);
//                System.out.println(i + ": " + startIdx + " " + best);
                int reach = 2 + i + possibleMixPieceSegments.get(startIdx + i + 2).size();
                if (reach > best) {
                    best = reach;
                    index = i;
                    bestAlpha = alphaTotal;
                } else if (reach == best) {
                    if (alphaTotal > bestAlpha) {
                        index = i;
                        bestAlpha = alphaTotal;
                    }
                }
            }
        }
        this.segments.add(possibleMixPieceSegments.get(startIdx).get(index));
        return startIdx + index + 1 + 1;
    }

    private List<MixPieceSegment> mergePerB(List<MixPieceSegment> segments) {
        double aMinTemp = -Double.MAX_VALUE;
        double aMaxTemp = Double.MAX_VALUE;
        double b = Double.NaN;
        List<Long> timestamps = new ArrayList<>();
        List<MixPieceSegment> mergedMixPieceSegments = new ArrayList<>();

        segments.sort(Comparator.comparingDouble(MixPieceSegment::getB).thenComparingDouble(MixPieceSegment::getA));
        for (int i = 0; i < segments.size(); i++) {
            if (b != segments.get(i).getB()) {
                if (timestamps.size() == 1)
                    mergedMixPieceSegments.add(new MixPieceSegment(timestamps.get(0), aMinTemp, aMaxTemp, b));
                else {
                    for (Long timestamp : timestamps)
                        mergedMixPieceSegments.add(new MixPieceSegment(timestamp, aMinTemp, aMaxTemp, b));
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
                if (timestamps.size() == 1) mergedMixPieceSegments.add(segments.get(i - 1));
                else {
                    for (long timestamp : timestamps)
                        mergedMixPieceSegments.add(new MixPieceSegment(timestamp, aMinTemp, aMaxTemp, b));
                }
                timestamps.clear();
                timestamps.add(segments.get(i).getInitTimestamp());
                aMinTemp = segments.get(i).getAMin();
                aMaxTemp = segments.get(i).getAMax();
            }
        }
        if (!timestamps.isEmpty()) {
            if (timestamps.size() == 1)
                mergedMixPieceSegments.add(new MixPieceSegment(timestamps.get(0), aMinTemp, aMaxTemp, b));
            else {
                for (long timestamp : timestamps)
                    mergedMixPieceSegments.add(new MixPieceSegment(timestamp, aMinTemp, aMaxTemp, b));
            }
        }

        return mergedMixPieceSegments;
    }


    private static void mergePerBNew(List<MixPieceSegment> segments, List<MixPieceSegment> mergedMixPieceSegments, List<MixPieceSegment> unmergedMixPieceSegments) {
        double aMinTemp = -Double.MAX_VALUE;
        double aMaxTemp = Double.MAX_VALUE;
        double b = Double.NaN;
        List<Long> timestamps = new ArrayList<>();

        segments.sort(Comparator.comparingDouble(MixPieceSegment::getB).thenComparingDouble(MixPieceSegment::getA));
        for (int i = 0; i < segments.size(); i++) {
            if (b != segments.get(i).getB()) {
                if (timestamps.size() == 1)
                    unmergedMixPieceSegments.add(new MixPieceSegment(timestamps.get(0), aMinTemp, aMaxTemp, b));
                else {
                    for (Long timestamp : timestamps)
                        mergedMixPieceSegments.add(new MixPieceSegment(timestamp, aMinTemp, aMaxTemp, b));
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
                if (timestamps.size() == 1) unmergedMixPieceSegments.add(segments.get(i - 1));
                else {
                    for (long timestamp : timestamps)
                        mergedMixPieceSegments.add(new MixPieceSegment(timestamp, aMinTemp, aMaxTemp, b));
                }
                timestamps.clear();
                timestamps.add(segments.get(i).getInitTimestamp());
                aMinTemp = segments.get(i).getAMin();
                aMaxTemp = segments.get(i).getAMax();
            }
        }
        if (!timestamps.isEmpty()) {
            if (timestamps.size() == 1)
                unmergedMixPieceSegments.add(new MixPieceSegment(timestamps.get(0), aMinTemp, aMaxTemp, b));
            else {
                for (long timestamp : timestamps)
                    mergedMixPieceSegments.add(new MixPieceSegment(timestamp, aMinTemp, aMaxTemp, b));
            }
        }
    }

    private static void mergeAll(List<MixPieceSegment> segments, List<MixPieceSegment> mergedMixPieceSegments, List<MixPieceSegment> unmergedMixPieceSegments) {
        double aMinTemp = -Double.MAX_VALUE;
        double aMaxTemp = Double.MAX_VALUE;
        List<Double> bValues = new ArrayList<>();
        List<Long> timestamps = new ArrayList<>();

        segments.sort(Comparator.comparingDouble(MixPieceSegment::getAMin));
        for (int i = 0; i < segments.size(); i++) {
            if (segments.get(i).getAMin() <= aMaxTemp && segments.get(i).getAMax() >= aMinTemp) {
                timestamps.add(segments.get(i).getInitTimestamp());
                aMinTemp = Math.max(aMinTemp, segments.get(i).getAMin());
                aMaxTemp = Math.min(aMaxTemp, segments.get(i).getAMax());
                bValues.add(segments.get(i).getB());
            } else {
                if (timestamps.size() == 1) unmergedMixPieceSegments.add(segments.get(i - 1));
                else {
                    for (int j = 0; j < timestamps.size(); j++)
                        mergedMixPieceSegments.add(new MixPieceSegment(timestamps.get(j), aMinTemp, aMaxTemp, bValues.get(j)));
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
                unmergedMixPieceSegments.add(new MixPieceSegment(timestamps.get(0), aMinTemp, aMaxTemp, bValues.get(0)));
            else {
                for (int i = 0; i < timestamps.size(); i++)
                    mergedMixPieceSegments.add(new MixPieceSegment(timestamps.get(i), aMinTemp, aMaxTemp, bValues.get(i)));
            }
        }
    }

    private static void merge(List<MixPieceSegment> segments, List<MixPieceSegment> perBMixPieceSegments, List<MixPieceSegment> perAMixPieceSegments, List<MixPieceSegment> restMixPieceSegments) {
        perBMixPieceSegments = new ArrayList<>();
        perAMixPieceSegments = new ArrayList<>();
        restMixPieceSegments = new ArrayList<>();
        List<MixPieceSegment> temp = new ArrayList<>();

        mergePerBNew(segments, perBMixPieceSegments, temp);
        if (!temp.isEmpty()) {
            mergeAll(temp, perAMixPieceSegments, restMixPieceSegments);
        }
    }

    public List<Point> decompress() {
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

    private void toByteArrayPerBSegments(List<MixPieceSegment> segments, boolean variableByte, ByteArrayOutputStream outStream) throws IOException {
        TreeMap<Integer, HashMap<Double, List<Long>>> input = new TreeMap<>();
        for (MixPieceSegment segment : segments) {
            double a = segment.getA();
            int b = (int) Math.round(segment.getB() / epsilon);
            long t = segment.getInitTimestamp();
            if (!input.containsKey(b)) input.put(b, new HashMap<>());
            if (!input.get(b).containsKey(a)) input.get(b).put(a, new ArrayList<>());
            input.get(b).get(a).add(t);
        }

        VariableByteEncoder.write(input.size(), outStream);
        if (input.isEmpty()) return;
        int previousB = input.firstKey();
        VariableByteEncoder.write(previousB, outStream);
        for (Map.Entry<Integer, HashMap<Double, List<Long>>> bMixPieceSegments : input.entrySet()) {
            VariableByteEncoder.write(bMixPieceSegments.getKey() - previousB, outStream);
            previousB = bMixPieceSegments.getKey();
            VariableByteEncoder.write(bMixPieceSegments.getValue().size(), outStream);
            for (Map.Entry<Double, List<Long>> aMixPieceSegment : bMixPieceSegments.getValue().entrySet()) {
                FloatEncoder.write(aMixPieceSegment.getKey().floatValue(), outStream);
                if (variableByte) Collections.sort(aMixPieceSegment.getValue());
                VariableByteEncoder.write(aMixPieceSegment.getValue().size(), outStream);
                long previousTS = 0;
                for (Long timestamp : aMixPieceSegment.getValue()) {
                    if (variableByte) VariableByteEncoder.write((int) (timestamp - previousTS), outStream);
                    else UIntEncoder.write(timestamp, outStream);
                    previousTS = timestamp;
                }
            }
        }
    }

    public byte[] toByteArray(boolean variableByte, boolean zstd) throws IOException {
        ByteArrayOutputStream outStream = new ByteArrayOutputStream();
        byte[] bytes;

        FloatEncoder.write((float) epsilon, outStream);
        toByteArrayPerBSegments(segments, variableByte, outStream);

        if (variableByte) VariableByteEncoder.write((int) lastTimeStamp, outStream);
        else UIntEncoder.write(lastTimeStamp, outStream);

        if (zstd) bytes = Zstd.compress(outStream.toByteArray());
        else bytes = outStream.toByteArray();

        outStream.close();

        return bytes;
    }

    private List<MixPieceSegment> readMergedPerBMixPieceSegments(boolean variableByte, ByteArrayInputStream inStream) throws IOException {
        ArrayList<MixPieceSegment> segments = new ArrayList<>();
        long numB = VariableByteEncoder.read(inStream);
        if (numB == 0) return segments;
        int previousB = VariableByteEncoder.read(inStream);
        for (int i = 0; i < numB; i++) {
            int b = VariableByteEncoder.read(inStream) + previousB;
            previousB = b;
            int numA = VariableByteEncoder.read(inStream);
            for (int j = 0; j < numA; j++) {
                float a = FloatEncoder.read(inStream);
                int numTimestamps = VariableByteEncoder.read(inStream);
                long timestamp = 0;
                for (int k = 0; k < numTimestamps; k++) {
                    if (variableByte) timestamp += VariableByteEncoder.read(inStream);
                    else timestamp = UIntEncoder.read(inStream);
                    segments.add(new MixPieceSegment(timestamp, a, (float) (b * epsilon)));
                }
            }
        }

        return segments;
    }

    private void readByteArray(byte[] input, boolean variableByte, boolean zstd) throws IOException {
        byte[] binary;
        if (zstd) binary = Zstd.decompress(input, input.length * 2); //TODO: How to know apriori original size?
        else binary = input;
        ByteArrayInputStream inStream = new ByteArrayInputStream(binary);

        this.epsilon = FloatEncoder.read(inStream);
        this.segments = readMergedPerBMixPieceSegments(variableByte, inStream);
        if (variableByte) this.lastTimeStamp = VariableByteEncoder.read(inStream);
        else this.lastTimeStamp = UIntEncoder.read(inStream);
        inStream.close();
    }
}
