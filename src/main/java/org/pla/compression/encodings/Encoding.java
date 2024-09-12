package org.pla.compression.encodings;

import com.github.benmanes.caffeine.cache.Caffeine;
import com.github.benmanes.caffeine.cache.LoadingCache;
import com.github.luben.zstd.Zstd;
import org.pla.compression.util.Encoding.FloatEncoder;
import org.pla.compression.util.Encoding.UIntEncoder;
import org.pla.compression.util.Encoding.VariableByteEncoder;
import org.pla.compression.util.Point;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.*;

public class Encoding {
    private static final double POWER = 0.1;
    private List<Segment> segments;

    private List<Segment> bestSegments;

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

    private double quantization(double value) {
        return Math.round(value / epsilon) * epsilon;
    }


    private double ceil(double value) {
        return Math.ceil(value / epsilon) * epsilon;
    }
    private double floor(double value) {
        return Math.floor(value / epsilon) * epsilon;
    }

    public List<Segment> getSegments() {
        return segments;
    }

    private List<Segment> createSegmentsFromStartIdx(int startIdx, List<Point> points) {
        List<Segment> segments = new LinkedList<>();

        long initTimestamp = points.get(startIdx).getTimestamp();
        double b1 = floor(points.get(startIdx).getValue());
        double b2 = ceil(points.get(startIdx).getValue());
        double b = quantization(points.get(startIdx).getValue());
        boolean floor = true;
        boolean ceil = true;
        int length = 0;
        if (startIdx + 1 == points.size()) {
            segments.add(new Segment(initTimestamp, -Double.MAX_VALUE, Double.MAX_VALUE, b));
            return segments;
        }
        double aMax1 = ((points.get(startIdx + 1).getValue() + epsilon) - b1) / (points.get(startIdx + 1).getTimestamp() - initTimestamp);
        double aMin1 = ((points.get(startIdx + 1).getValue() - epsilon) - b1) / (points.get(startIdx + 1).getTimestamp() - initTimestamp);
        double aMax2 = ((points.get(startIdx + 1).getValue() + epsilon) - b2) / (points.get(startIdx + 1).getTimestamp() - initTimestamp);
        double aMin2 = ((points.get(startIdx + 1).getValue() - epsilon) - b2) / (points.get(startIdx + 1).getTimestamp() - initTimestamp);

        if (startIdx + 2 == points.size()) {
            segments.add(new Segment(initTimestamp, aMin1, aMax1, b1));
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
                segments.add(new Segment(initTimestamp, aMin1, aMax1, b1));
            } else if (length < 0){
                segments.add(new Segment(initTimestamp, aMin2, aMax2, b2));
            } else if (aMax1 - aMin1 > aMax2 - aMin2) {
                segments.add(new Segment(initTimestamp, aMin1, aMax1, b1));
            } else {
                segments.add(new Segment(initTimestamp, aMin2, aMax2, b2));
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
            segments.add(new Segment(initTimestamp, aMin1, aMax1, b1));
        } else if (length < 0){
            segments.add(new Segment(initTimestamp, aMin2, aMax2, b2));
        } else if (aMax1 - aMin1 > aMax2 - aMin2) {
            segments.add(new Segment(initTimestamp, aMin1, aMax1, b1));
        } else {
            segments.add(new Segment(initTimestamp, aMin2, aMax2, b2));
        }

        return segments;
    }

    private int findReachFrom(int startIdx, List<Point> points) {
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
        double b1 = floor(points.get(startIdx).getValue());
        double b2 = ceil(points.get(startIdx).getValue());
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


    private State binarySearch(List<State> states, int previousNoOfGroups, long initTimestamp, double b, List<Segment> segments) {
        List<Segment> oldSegments = mergePerB(segments);
        TreeMap<Double, HashMap<Double, List<Long>>> tree = new TreeMap<>();
        for (Segment segment : oldSegments) {
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
            if (findNumberOfGroups(initTimestamp, state.aMax, state.aMin, b, segments, previousNoOfGroups, oldSegments, tree)) {
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
    private int findNumberOfGroups(long initTimestamp, double aMax, double aMin, double b, List<Segment> segments) {
        List<Segment> tempSegments = new ArrayList<>(segments);
        tempSegments.add(new Segment(initTimestamp, aMin, aMax, b));
        tempSegments = mergePerB(tempSegments);

        TreeMap<Double, HashMap<Double, ArrayList<Long>>> input = new TreeMap<>();
        for (Segment segment : tempSegments) {
            double a = segment.getA();
            double bi = segment.getB();
            long t = segment.getInitTimestamp();
            if (!input.containsKey(bi)) input.put(bi, new HashMap<>());
            if (!input.get(bi).containsKey(a)) input.get(bi).put(a, new ArrayList<>());
            input.get(bi).get(a).add(t);
        }
        int groups = 0;
        for (Map.Entry<Double, HashMap<Double, ArrayList<Long>>> bSegments : input.entrySet())
            groups += bSegments.getValue().size();
        return groups;
    }

    private boolean findNumberOfGroups(long initTimestamp, double aMax, double aMin, double b, List<Segment> segments,
                                       int previousNoOfGroups, List<Segment> oldSegments,
                                       TreeMap<Double, HashMap<Double, List<Long>>> tree) {
        double alpha = oldSegments.stream().filter(s -> s.getB() == b && s.getAMin() <= aMax && s.getAMax() >= aMin).mapToDouble(s -> s.getAMax() - s.getAMin()).max().orElse(-1.0);

        List<Segment> tempSegments = new ArrayList<>(segments);
        tempSegments.add(new Segment(initTimestamp, aMin, aMax, b));
        tempSegments = mergePerB(tempSegments);
        double value = 0;
        TreeMap<Double, HashMap<Double, List<Long>>> input = new TreeMap<>();
        for (Segment segment : tempSegments) {
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
        for (Map.Entry<Double, HashMap<Double, List<Long>>> bSegments : input.entrySet())
            groups += bSegments.getValue().size();

        return groups > previousNoOfGroups;// || value > 10;
    }

    LoadingCache<Integer, List<Segment>> cache = Caffeine.newBuilder()
            .maximumSize(5000)
            .build(this::createSegmentsFromStartIdx);
//    Map<Integer, List<SimPieceSegment>> cache = new HashMap<>();

    private List<Segment> createSegmentsFromStartIdx(int idx) {
        return createSegmentsFromStartIdx(idx, this.points);
    }

    private List<Segment> compress(List<Point> points) {
        return compress(points, 0, 0);
    }

    private List<Segment> compress(List<Point> points, int mode, double pow) {
        Map<Integer, List<Segment>> possibleSegments = new TreeMap<>();
        switch (mode) {
            case 0:
                for (int i = 0; i <= points.size() - 1; i++) {
                    List<Segment> segmentsFromStartIdx = createSegmentsFromStartIdx(i, points);
//            System.out.println(String.format("StartIdx: %d: %d", i, segmentsFromStartIdx.size()));
                    possibleSegments.put(i, segmentsFromStartIdx);
                }
                this.segments = new ArrayList<>();
                double[][] best = new double[points.size()][];
                double angle = possibleSegments.get(points.size() - 1).get(0).getAMax() - possibleSegments.get(points.size() - 1).get(0).getAMax();
                best[points.size() - 1] = new double[]{1, Math.pow(angle, pow), 1};
                for (int i=points.size()-2; i>=0; i--) {
                    findBestWithAngle(i, possibleSegments, best, pow);
                }

                int start = 0;
                int count = 0;
                while (start < points.size()) {
                    this.segments.add(possibleSegments.get(start).get((int) (best[start][2]-1)));
                    start += (int) (best[start][2]) + 1;
                    count++;
                }

//                System.out.println("Best: " + best[0][0] + " Used: " + this.segments.size());

                this.segments = mergePerB(this.segments);
                break;
            case 1:
                List<Segment> segments = new ArrayList<>();
                int currentIdx = 0;
                while (currentIdx < points.size()) currentIdx = createSegment(currentIdx, points, segments);
                this.segments = mergePerB(segments);
                break;
            case 2:
                this.segments = new ArrayList<>();
                int startIdx = 0;
                while (startIdx < points.size()) {
                    startIdx = addSegment(startIdx, pow);
                }
                this.segments = mergePerB(this.segments);
                break;
            case 3:
                for (int i = 0; i <= points.size() - 1; i++) {
                    List<Segment> segmentsFromStartIdx = createSegmentsFromStartIdx(i, points);
                    possibleSegments.put(i, segmentsFromStartIdx);
                }
                this.segments = new ArrayList<>();
                int startIdxG = 0;
                while (startIdxG < points.size()) {
                    startIdxG = addSegment(startIdxG, possibleSegments, pow);
                }
                this.segments = mergePerB(this.segments);
                break;
            case 4:
                for (int i = 0; i <= points.size() - 1; i++) {
                    List<Segment> segmentsFromStartIdx = createSegmentsFromStartIdx(i, points);
                    possibleSegments.put(i, segmentsFromStartIdx);
                }
                this.segments = new ArrayList<>();
                startIdx = 0;
                while (startIdx < points.size()) {
                    startIdx = addSegment2(startIdx, pow);
                }
                this.segments = mergePerB(this.segments);
                break;
            case 5:
                List<Segment> segmentsBoth = new ArrayList<>();
                int currentIdxBoth = 0;
                while (currentIdxBoth < points.size()) currentIdxBoth = createSegmentBoth(currentIdxBoth, points, segmentsBoth);
                this.segments = mergePerB(segmentsBoth);
                break;
        }
        return segments;
    }
    
     private byte[] compressMixPiece(List<Point> points, int mode, double pow) {
        Map<Integer, List<Segment>> possibleSegments = new TreeMap<>();
        switch (mode) {
            case 5:
                List<Segment> segmentsBoth = new ArrayList<>();
                int currentIdxBoth = 0;
                while (currentIdxBoth < points.size()) currentIdxBoth = createSegmentBoth(currentIdxBoth, points, segmentsBoth);
                int globalMinB = (int)(segments.stream().mapToDouble(Segment::getB).min().orElse(Integer.MIN_VALUE) / epsilon);
                List<Segment> perBSegments = new ArrayList<Segment>();
                List<Segment> perASegments = new ArrayList<Segment>();
                List<Segment> restSegments = new ArrayList<Segment>();
                merge(segmentsBoth, perBSegments, perASegments, restSegments);
                return toByteArray(epsilon, globalMinB, perBSegments, perASegments, restSegments);
        }
        return null;
    }
    
    private byte[] toByteArray(double epsilon, int globalMinB, List<Segment> perBSegments, List<Segment> perASegments, List<Segment> restSegments) {
        ByteArrayOutputStream outStream = new ByteArrayOutputStream();
        byte[] bytes = null;

        try {
            FloatEncoder.write((float) epsilon, outStream);
            VariableByteEncoder.write(globalMinB, outStream);

            toByteArrayPerBSegments(perBSegments, outStream, epsilon, globalMinB);
            toByteArrayPerASegments(perASegments, outStream, epsilon, globalMinB);
            toByteArrayRestSegments(restSegments, outStream, epsilon, globalMinB);

            VariableByteEncoder.write((int) lastTimeStamp, outStream);

             return outStream.toByteArray();
            //bytes = Zstd.compress(outStream.toByteArray());

            //outStream.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        return bytes;
    }


    private int createSegment(int startIdx, List<Point> points, List<Segment> segments) {
        long initTimestamp = points.get(startIdx).getTimestamp();
        double b = floor(points.get(startIdx).getValue());
        if (startIdx + 1 == points.size()) {
            segments.add(new Segment(initTimestamp, -Double.MAX_VALUE, Double.MAX_VALUE, b));
            return startIdx + 1;
        }
        double aMax = ((points.get(startIdx + 1).getValue() + epsilon) - b) / (points.get(startIdx + 1).getTimestamp() - initTimestamp);
        double aMin = ((points.get(startIdx + 1).getValue() - epsilon) - b) / (points.get(startIdx + 1).getTimestamp() - initTimestamp);
        if (startIdx + 2 == points.size()) {
            segments.add(new Segment(initTimestamp, aMin, aMax, b));
            return startIdx + 2;
        }

        for (int idx = startIdx + 2; idx < points.size(); idx++) {
            double upValue = points.get(idx).getValue() + epsilon;
            double downValue = points.get(idx).getValue() - epsilon;

            double upLim = aMax * (points.get(idx).getTimestamp() - initTimestamp) + b;
            double downLim = aMin * (points.get(idx).getTimestamp() - initTimestamp) + b;
            if ((downValue > upLim || upValue < downLim)) {
                segments.add(new Segment(initTimestamp, aMin, aMax, b));
                return idx;
            }

            if (upValue < upLim)
                aMax = Math.max((upValue - b) / (points.get(idx).getTimestamp() - initTimestamp), aMin);
            if (downValue > downLim)
                aMin = Math.min((downValue - b) / (points.get(idx).getTimestamp() - initTimestamp), aMax);
        }
        segments.add(new Segment(initTimestamp, aMin, aMax, b));

        return points.size();
    }

    private int createSegmentBoth(int startIdx, List<Point> points, List<Segment> segments) {
        long initTimestamp = points.get(startIdx).getTimestamp();
        double b1 = floor(points.get(startIdx).getValue());
        double b2 = ceil(points.get(startIdx).getValue());
        int count = 0;
        boolean floor = true;
        boolean ceil = true;
        double b = quantization(points.get(startIdx).getValue());
        if (startIdx + 1 == points.size()) {
            segments.add(new Segment(initTimestamp, -Double.MAX_VALUE, Double.MAX_VALUE, b));
            return startIdx + 1;
        }
        double aMax1 = ((points.get(startIdx + 1).getValue() + epsilon) - b1) / (points.get(startIdx + 1).getTimestamp() - initTimestamp);
        double aMin1 = ((points.get(startIdx + 1).getValue() - epsilon) - b1) / (points.get(startIdx + 1).getTimestamp() - initTimestamp);
        double aMax2 = ((points.get(startIdx + 1).getValue() + epsilon) - b2) / (points.get(startIdx + 1).getTimestamp() - initTimestamp);
        double aMin2 = ((points.get(startIdx + 1).getValue() - epsilon) - b2) / (points.get(startIdx + 1).getTimestamp() - initTimestamp);
        if (startIdx + 2 == points.size()) {
            segments.add(new Segment(initTimestamp, aMin1, aMax1, b));
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
                    segments.add(new Segment(initTimestamp, aMin1, aMax1, b1));
                } else {
                    segments.add(new Segment(initTimestamp, aMin2, aMax2, b2));
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
            segments.add(new Segment(initTimestamp, aMin1, aMax1, b1));
        } else {
            segments.add(new Segment(initTimestamp, aMin2, aMax2, b2));
        }

        return points.size();
    }



    private int createSegmentSimPath(int startIdx, List<Point> points, List<Segment> segments) {
//        System.out.println("Starting segment " + segments.size());
        long initTimestamp = points.get(startIdx).getTimestamp();
        double b = quantization(points.get(startIdx).getValue());
        if (startIdx + 1 == points.size()) {
            segments.add(new Segment(initTimestamp, -Double.MAX_VALUE, Double.MAX_VALUE, b));
            return startIdx + 1;
        }
        double aMax = ((points.get(startIdx + 1).getValue() + epsilon) - b) / (points.get(startIdx + 1).getTimestamp() - initTimestamp);
        double aMin = ((points.get(startIdx + 1).getValue() - epsilon) - b) / (points.get(startIdx + 1).getTimestamp() - initTimestamp);
//      System.out.println(aMax-aMin);
        double diff = aMax-aMin;
        if (startIdx + 2 == points.size()) {
            segments.add(new Segment(initTimestamp, aMin, aMax, b));
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
//                    System.out.println(String.format("Segments: %d, Groups: %d, Size: %d", segments.size(), newGroups, (idx - startIdx)));
                    State state = binarySearch(states, previousNoOfGroups, initTimestamp, b, segments);
//                    if (idx > state.idx) {
//                        System.out.println(String.format("%d - %d - %f", (idx - startIdx), (state.idx - startIdx) , ((double) state.idx - startIdx) / (idx - startIdx)));
//                    }
                    previousNoOfGroups = newGroups;
                    segments.add(new Segment(initTimestamp, state.aMin, state.aMax, b));
                    return state.idx;
                }
                previousNoOfGroups = newGroups;
                segments.add(new Segment(initTimestamp, aMin, aMax, b));
                return idx;
            }
            diff = aMax-aMin;
            if (upValue < upLim)
                aMax = Math.max((upValue - b) / (points.get(idx).getTimestamp() - initTimestamp), aMin);
            if (downValue > downLim)
                aMin = Math.min((downValue - b) / (points.get(idx).getTimestamp() - initTimestamp), aMax);
            diff = diff - (aMax-aMin);
        }
        segments.add(new Segment(initTimestamp, aMin, aMax, b));

        return points.size();
    }

    private int findBest(int start, Map<Integer, List<Segment>> possibleSegments, int[][] best) {
//        System.out.println("CALLED: " + start + " - " + best.length);
        if (start >= possibleSegments.size()) {
            return 0;
        }
        if (best[start] != null) {
            return best[start][0];
        }
        else {
            int bestResult = Integer.MAX_VALUE;
            int bestIndex = 0;
            for (int i=1; i<= possibleSegments.get(start).size(); i++) {
                int result = 1 + findBest(start + i + 1, possibleSegments, best);
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



    private double[] findBestWithAngle(int start, Map<Integer, List<Segment>> possibleSegments, double[][] best, double pow) {
        if (start >= possibleSegments.size()) {
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
            for (int i=1; i<= possibleSegments.get(start).size(); i++) {
                double[] result = findBestWithAngle(start + i + 1, possibleSegments, best, pow);
                double n = result[0] + 1;
//                System.out.println(start);
                Segment segment = possibleSegments.get(start).get(i-1);
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
//            System.out.println("Best for: " + start + " is: " + bestResult);
            return new double[] {bestN, bestAngle};
        }
    }

    private int addSegment(int startIdx, double pow) {
//        int firstSegments = cache.get(startIdx).size();
        List<Segment> firstSegments = createSegmentsFromStartIdx(startIdx);
        int index = 0;
        if (startIdx + 2 < points.size()) {
//            int best = 2 + cache.get(startIdx + 2).size();
            int best = 2 + findReachFrom(startIdx + 2, points);
            for (int i = (int) Math.pow(firstSegments.size(), pow); i<firstSegments.size() && startIdx + i + 2 < points.size(); i++) {
//                int reach = 2 + i + cache.get(startIdx + i + 2).size();
                int reach = 2 + i + findReachFrom(startIdx + i + 2, points);
                if (reach > best) {
                    best = reach;
                    index = i;
                }
            }
        }
//        this.segments.add(cache.get(startIdx).get(index));
        this.segments.add(firstSegments.get(index));
        return startIdx + index + 1 + 1;
    }


    private int addSegment2(int startIdx, double pow) {


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

        List<Segment> firstSegments = createSegmentsFromStartIdx(startIdx);
        Segment currSegment = firstSegments.get(firstSegments.size() - 1);
        if (!this.segments.isEmpty()) {
            Segment prevSegment = this.segments.get(this.segments.size() - 1);
            List<Segment> startPrevSegments = createSegmentsFromStartIdx((int) prevSegment.getInitTimestamp());
            int bLast = startIdx + firstSegments.size();
            int bPrev = startIdx - 1;
            int bStart = (int) prevSegment.getInitTimestamp();
            double sumAngles = Math.pow(prevSegment.getAMax() - prevSegment.getAMin(), pow) + Math.pow(currSegment.getAMax() - currSegment.getAMin(), pow);
            int bestPrev = bPrev;
            while (bPrev >= 0) {
                List<Segment> prevSegments = createSegmentsFromStartIdx(bPrev);
                if (bPrev + prevSegments.size() + 1 < bLast) {
                    break;
                }
                if (bLast - bPrev - 1 >= prevSegments.size() || bPrev - bStart - 1 - 1 < 1) {
                    bPrev--;
                    continue;
                }
                Segment newLastSegment = prevSegments.get(bLast - bPrev - 1);
                Segment newPrevSegment = startPrevSegments.get(bPrev - bStart - 1 - 1);
                double newSum = Math.pow(newLastSegment.getAMax() - newLastSegment.getAMin(), pow) + Math.pow(newPrevSegment.getAMax() - newPrevSegment.getAMin(), pow);
                if (newSum > sumAngles) {
                    sumAngles = newSum;
                    bestPrev = bPrev;
                    prevSegment = newPrevSegment;
                    currSegment = newLastSegment;
                }
                bPrev--;
            }
            this.segments.remove(this.segments.size()-1);
            this.segments.add(prevSegment);
        }
        this.segments.add(currSegment);
        return startIdx + firstSegments.size() + 1;
    }


    private int addSegment(int startIdx, Map<Integer, List<Segment>> possibleSegments, double pow) {
//        System.out.println("Getting: " + startIdx);
//        System.out.println("Size: " + possibleSegments.get(startIdx).size());
        int firstSegments = possibleSegments.get(startIdx).size();
        int index = 0;
        if (startIdx + 2 < possibleSegments.size()) {
            double alpha1 = possibleSegments.get(startIdx).get(0).getAMax() - possibleSegments.get(startIdx).get(0).getAMax();
            double bestAlpha = Math.pow(alpha1, pow) + Math.pow(
                    possibleSegments.get(startIdx + 2).get(possibleSegments.get(startIdx + 2).size() -1).getAMax() -
                            possibleSegments.get(startIdx + 2).get(possibleSegments.get(startIdx + 2).size() -1).getAMin(), pow);
            int best = 2 + possibleSegments.get(startIdx + 2).size();
//            System.out.println(0 + ": " + startIdx + " " + best);
            for (int i=1; i<firstSegments && startIdx + i + 2 < possibleSegments.size(); i++) {
                double alphaTotal = Math.pow(alpha1, pow) + Math.pow(
                        possibleSegments.get(startIdx + i + 2).get(possibleSegments.get(startIdx + i + 2).size() -1).getAMax() -
                                possibleSegments.get(startIdx + i + 2).get(possibleSegments.get(startIdx + i + 2).size() -1).getAMin(), pow);
//                System.out.println(i + ": " + startIdx + " " + best);
                int reach = 2 + i + possibleSegments.get(startIdx + i + 2).size();
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
        this.segments.add(possibleSegments.get(startIdx).get(index));
        return startIdx + index + 1 + 1;
    }
    
    private List<Segment> mergePerB(List<Segment> segments) {
        double aMinTemp = -Double.MAX_VALUE;
        double aMaxTemp = Double.MAX_VALUE;
        double b = Double.NaN;
        List<Long> timestamps = new ArrayList<>();
        List<Segment> mergedSegments = new ArrayList<>();

        segments.sort(Comparator.comparingDouble(Segment::getB).thenComparingDouble(Segment::getA));
        for (int i = 0; i < segments.size(); i++) {
            if (b != segments.get(i).getB()) {
                if (timestamps.size() == 1)
                    mergedSegments.add(new Segment(timestamps.get(0), aMinTemp, aMaxTemp, b));
                else {
                    for (Long timestamp : timestamps)
                        mergedSegments.add(new Segment(timestamp, aMinTemp, aMaxTemp, b));
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
                if (timestamps.size() == 1) mergedSegments.add(segments.get(i - 1));
                else {
                    for (long timestamp : timestamps)
                        mergedSegments.add(new Segment(timestamp, aMinTemp, aMaxTemp, b));
                }
                timestamps.clear();
                timestamps.add(segments.get(i).getInitTimestamp());
                aMinTemp = segments.get(i).getAMin();
                aMaxTemp = segments.get(i).getAMax();
            }
        }
        if (!timestamps.isEmpty()) {
            if (timestamps.size() == 1)
                mergedSegments.add(new Segment(timestamps.get(0), aMinTemp, aMaxTemp, b));
            else {
                for (long timestamp : timestamps)
                    mergedSegments.add(new Segment(timestamp, aMinTemp, aMaxTemp, b));
            }
        }

        return mergedSegments;
    }
    
    
    private static void mergePerBNew(List<Segment> segments, List<Segment> mergedSegments, List<Segment> unmergedSegments) {
        double aMinTemp = -Double.MAX_VALUE;
        double aMaxTemp = Double.MAX_VALUE;
        double b = Double.NaN;
        List<Long> timestamps = new ArrayList<>();

        segments.sort(Comparator.comparingDouble(Segment::getB).thenComparingDouble(Segment::getA));
        for (int i = 0; i < segments.size(); i++) {
            if (b != segments.get(i).getB()) {
                if (timestamps.size() == 1)
                    unmergedSegments.add(new Segment(timestamps.get(0), aMinTemp, aMaxTemp, b));
                else {
                    for (Long timestamp : timestamps)
                        mergedSegments.add(new Segment(timestamp, aMinTemp, aMaxTemp, b));
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
                        mergedSegments.add(new Segment(timestamp, aMinTemp, aMaxTemp, b));
                }
                timestamps.clear();
                timestamps.add(segments.get(i).getInitTimestamp());
                aMinTemp = segments.get(i).getAMin();
                aMaxTemp = segments.get(i).getAMax();
            }
        }
        if (!timestamps.isEmpty()) {
            if (timestamps.size() == 1)
                unmergedSegments.add(new Segment(timestamps.get(0), aMinTemp, aMaxTemp, b));
            else {
                for (long timestamp : timestamps)
                    mergedSegments.add(new Segment(timestamp, aMinTemp, aMaxTemp, b));
            }
        }
    }

    private static void mergeAll(List<Segment> segments, List<Segment> mergedSegments, List<Segment> unmergedSegments) {
        double aMinTemp = -Double.MAX_VALUE;
        double aMaxTemp = Double.MAX_VALUE;
        List<Double> bValues = new ArrayList<>();
        List<Long> timestamps = new ArrayList<>();

        segments.sort(Comparator.comparingDouble(Segment::getAMin));
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
                        mergedSegments.add(new Segment(timestamps.get(j), aMinTemp, aMaxTemp, bValues.get(j)));
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
                unmergedSegments.add(new Segment(timestamps.get(0), aMinTemp, aMaxTemp, bValues.get(0)));
            else {
                for (int i = 0; i < timestamps.size(); i++)
                    mergedSegments.add(new Segment(timestamps.get(i), aMinTemp, aMaxTemp, bValues.get(i)));
            }
        }
    }

    private static void merge(List<Segment> segments, List<Segment> perBSegments, List<Segment> perASegments, List<Segment> restSegments) {
        perBSegments = new ArrayList<>();
        perASegments = new ArrayList<>();
        restSegments = new ArrayList<>();
        List<Segment> temp = new ArrayList<>();

        mergePerBNew(segments, perBSegments, temp);
        if (!temp.isEmpty()) {
            mergeAll(temp, perASegments, restSegments);
        }
    }

    public List<Point> decompress() {
        segments.sort(Comparator.comparingLong(Segment::getInitTimestamp));
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

    private void toByteArrayPerBSegments(List<Segment> segments, boolean variableByte, ByteArrayOutputStream outStream) throws IOException {
        TreeMap<Integer, HashMap<Double, List<Long>>> input = new TreeMap<>();
        for (Segment segment : segments) {
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
        for (Map.Entry<Integer, HashMap<Double, List<Long>>> bSegments : input.entrySet()) {
            VariableByteEncoder.write(bSegments.getKey() - previousB, outStream);
            previousB = bSegments.getKey();
            VariableByteEncoder.write(bSegments.getValue().size(), outStream);
            for (Map.Entry<Double, List<Long>> aSegment : bSegments.getValue().entrySet()) {
                FloatEncoder.write(aSegment.getKey().floatValue(), outStream);
                if (variableByte) Collections.sort(aSegment.getValue());
                VariableByteEncoder.write(aSegment.getValue().size(), outStream);
                long previousTS = 0;
                for (Long timestamp : aSegment.getValue()) {
                    if (variableByte) VariableByteEncoder.write((int) (timestamp - previousTS), outStream);
                    else UIntEncoder.write(timestamp, outStream);
                    previousTS = timestamp;
                }
            }
        }
    }

    private static void toByteArrayPerBSegments(List<Segment> segments, ByteArrayOutputStream outStream, double epsilon, int globalMinB) throws IOException {
        TreeMap<Integer, HashMap<Double, ArrayList<Long>>> input = new TreeMap<>();
        for (Segment segment : segments) {
            double a = segment.getA();
            int b = (int) Math.round(segment.getB() / epsilon);
            long t = segment.getInitTimestamp();
            if (!input.containsKey(b)) input.put(b, new HashMap<>());
            if (!input.get(b).containsKey(a)) input.get(b).put(a, new ArrayList<>());
            input.get(b).get(a).add(t);
        }

        VariableByteEncoder.write(input.size(), outStream);
        if (input.isEmpty())
            return;
        int previousB = input.firstKey() - globalMinB;
        VariableByteEncoder.write(previousB, outStream);
        for (Map.Entry<Integer, HashMap<Double, ArrayList<Long>>> bSegments : input.entrySet()) {
            VariableByteEncoder.write(bSegments.getKey() - globalMinB - previousB, outStream);
            previousB = bSegments.getKey() - globalMinB;
            VariableByteEncoder.write(bSegments.getValue().size(), outStream);
            for (Map.Entry<Double, ArrayList<Long>> aSegment : bSegments.getValue().entrySet()) {
                FloatEncoder.write(aSegment.getKey().floatValue(), outStream);
                Collections.sort(aSegment.getValue());
                VariableByteEncoder.write(aSegment.getValue().size(), outStream);
                long previousTS = 0;
                for (Long timestamp : aSegment.getValue()) {
                    VariableByteEncoder.write((int) (timestamp - previousTS), outStream);
                    previousTS = timestamp;
                }
            }
        }
    }

    private static void toByteArrayPerASegments(List<Segment> segments, ByteArrayOutputStream outStream, double epsilon, int globalMinB) throws IOException {
        TreeMap<Double, List<Segment>> input = new TreeMap<>();
        for (Segment segment : segments) {
            if (!input.containsKey(segment.getA())) input.put(segment.getA(), new ArrayList<>());
            input.get(segment.getA()).add(segment);
        }

        VariableByteEncoder.write(input.size(), outStream);
        for (Map.Entry<Double, List<Segment>> aSegments : input.entrySet()) {
            FloatEncoder.write(aSegments.getKey().floatValue(), outStream);
            VariableByteEncoder.write(aSegments.getValue().size(), outStream);
            aSegments.getValue().sort(Comparator.comparingDouble(Segment::getB));
            int previousB = (int) Math.round(aSegments.getValue().get(0).getB() / epsilon) - globalMinB;
            VariableByteEncoder.write(previousB, outStream);
            for (Segment segment : aSegments.getValue()) {
                VariableByteEncoder.write((int) (Math.round(segment.getB() / epsilon) - globalMinB - previousB), outStream);
                previousB = (int) Math.round(segment.getB() / epsilon) - globalMinB;
                UIntEncoder.write(segment.getInitTimestamp(), outStream);
            }
        }
    }

    private static void toByteArrayRestSegments(List<Segment> segments, ByteArrayOutputStream outStream, double epsilon, int globalMinB) throws IOException {
        VariableByteEncoder.write(segments.size(), outStream);
        if (segments.isEmpty())
            return;
        segments.sort(Comparator.comparingDouble(Segment::getB));
        int previousB = (int) Math.round(segments.get(0).getB() / epsilon) - globalMinB;
        VariableByteEncoder.write(previousB, outStream);
        for (Segment segment : segments) {
            VariableByteEncoder.write((int) (Math.round(segment.getB() / epsilon) - globalMinB - previousB), outStream);
            previousB = (int) Math.round(segment.getB() / epsilon) - globalMinB;
            FloatEncoder.write((float) segment.getA(), outStream);
            UIntEncoder.write(segment.getInitTimestamp(), outStream);
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

    private List<Segment> readMergedPerBSegments(boolean variableByte, ByteArrayInputStream inStream) throws IOException {
        ArrayList<Segment> segments = new ArrayList<>();
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
                    segments.add(new Segment(timestamp, a, (float) (b * epsilon)));
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
        this.segments = readMergedPerBSegments(variableByte, inStream);
        if (variableByte) this.lastTimeStamp = VariableByteEncoder.read(inStream);
        else this.lastTimeStamp = UIntEncoder.read(inStream);
        inStream.close();
    }
}
