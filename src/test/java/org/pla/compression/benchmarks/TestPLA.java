package org.pla.compression.benchmarks;

import org.pla.compression.encodings.Encoding;
import org.pla.compression.encodings.MixPiece;
import org.pla.compression.encodings.Point;
import org.pla.compression.encodings.Result;
import org.pla.compression.util.TimeSeries;
import org.pla.compression.util.TimeSeriesReader;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.time.Duration;
import java.time.Instant;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class TestPLA {
    private Duration compressDuration;
    private Duration decompressDuration;

    private double[] Encoding(List<Point> ts, double epsilon, boolean variableByte, boolean zstd, int mode, double pow) throws IOException {
        compressDuration = Duration.ZERO;
        Instant start = Instant.now();
        Encoding encoding = new Encoding(ts, epsilon, mode, pow);
        compressDuration = Duration.between(start, Instant.now());
        byte[] binary = encoding.toByteArray(variableByte, zstd);

        encoding = new Encoding(binary, variableByte, zstd);
        Instant decompress = Instant.now();
        List<Point> tsDecompressed = encoding.decompress();
        decompressDuration = Duration.between(decompress, Instant.now());
        int idx = 0;
        double ae = 0.0;
        double mse = 0.0;
        for (Point expected : tsDecompressed) {
            Point actual = ts.get(idx);
            if (expected.getTimestamp() != actual.getTimestamp()) continue;
            idx++;
            double diff = actual.getValue() - expected.getValue();
            ae += Math.abs(diff);
            mse += Math.abs(diff*diff);
            assertEquals(actual.getValue(), expected.getValue(), 1.1 * epsilon, "Value did not match for timestamp " + actual.getTimestamp());
        }
        assertEquals(idx, ts.size());

        return new double[]{binary.length, encoding.getSegments().size(), (long) ae, mse};
    }

    private double[] MixPiece(List<Point> ts, double epsilon, boolean variableByte, boolean zstd) throws Exception {
        compressDuration = Duration.ZERO;
        Instant start = Instant.now();
        Result result = MixPiece.compress(ts, epsilon, variableByte, zstd);
        compressDuration = Duration.between(start, Instant.now());
        Instant decompress = Instant.now();
        List<Point> tsDecompressed = MixPiece.decompress(result.getBinary(), variableByte, zstd);
        decompressDuration = Duration.between(decompress, Instant.now());
        int idx = 0;
        double ae = 0.0;
        double mse = 0.0;
        for (Point expected : tsDecompressed) {
            Point actual = ts.get(idx);
            if (expected.getTimestamp() != actual.getTimestamp()) continue;
            idx++;
            double diff = actual.getValue() - expected.getValue();
            ae += Math.abs(diff);
            mse += Math.abs(diff*diff);
            assertEquals(actual.getValue(), expected.getValue(), 1.1 * epsilon, "Value did not match for timestamp " + actual.getTimestamp());
        }
        assertEquals(idx, ts.size());

        return new double[]{result.getBinary().length, result.getSegmentsSize(), (long) ae, mse};
    }

    private double[] MixPieceTunablePeekAhead(List<Point> ts, double epsilon, boolean variableByte, boolean zstd, double pow) throws Exception {
        compressDuration = Duration.ZERO;
        Instant start = Instant.now();
        Result result = MixPiece.compressTunablePeekAhead(ts, epsilon, variableByte, zstd, pow);
        compressDuration = Duration.between(start, Instant.now());
        Instant decompress = Instant.now();
        List<Point> tsDecompressed = null;
        for (int i=0;i<5;i++) {
            tsDecompressed = MixPiece.decompressImproved(result.getBinary(), variableByte, zstd);

        }
        decompressDuration = Duration.between(decompress, Instant.now()).dividedBy(5);
        int idx = 0;
        double ae = 0.0;
        double mse = 0.0;
        for (Point expected : tsDecompressed) {
            Point actual = ts.get(idx);
            if (expected.getTimestamp() != actual.getTimestamp()) continue;
            idx++;
            double diff = actual.getValue() - expected.getValue();
            ae += Math.abs(diff);
            mse += Math.abs(diff*diff);
            assertEquals(actual.getValue(), expected.getValue(), 1.1 * epsilon, "Value did not match for timestamp " + actual.getTimestamp());
        }
        assertEquals(idx, ts.size());

        return new double[]{result.getBinary().length, result.getSegmentsSize(), (long) ae, mse};
    }

    private double[] MixPieceQuantOptimal(List<Point> ts, double epsilon, boolean variableByte, boolean zstd, double pow) throws Exception {
        compressDuration = Duration.ZERO;
        Instant start = Instant.now();
        Result result = MixPiece.compressQuantOptimal(ts, epsilon, variableByte, zstd, pow);
        compressDuration = Duration.between(start, Instant.now());
        Instant decompress = Instant.now();
        List<Point> tsDecompressed = null;
        for (int i=0;i<5;i++) {
            tsDecompressed = MixPiece.decompressImproved(result.getBinary(), variableByte, zstd);

        }
        decompressDuration = Duration.between(decompress, Instant.now()).dividedBy(5);
        int idx = 0;
        double ae = 0.0;
        double mse = 0.0;
        for (Point expected : tsDecompressed) {
            Point actual = ts.get(idx);
            if (expected.getTimestamp() != actual.getTimestamp()) continue;
            idx++;
            double diff = actual.getValue() - expected.getValue();
            ae += Math.abs(diff);
            mse += Math.abs(diff*diff);
            assertEquals(actual.getValue(), expected.getValue(), 1.1 * epsilon, "Value did not match for timestamp " + actual.getTimestamp());
        }
        assertEquals(idx, ts.size());

        return new double[]{result.getBinary().length, result.getSegmentsSize(), (long) ae, mse};
    }


    private double[] MixPieceMinimumAngle(List<Point> ts, double epsilon, boolean variableByte, boolean zstd) throws Exception {
        compressDuration = Duration.ZERO;
        Instant start = Instant.now();
        Result result = MixPiece.compressMinimumAngle(ts, epsilon, variableByte, zstd);
        compressDuration = Duration.between(start, Instant.now());
        Instant decompress = Instant.now();
        List<Point> tsDecompressed = null;
        for (int i=0;i<5;i++) {
            tsDecompressed = MixPiece.decompressImproved(result.getBinary(), variableByte, zstd);

        }
        decompressDuration = Duration.between(decompress, Instant.now()).dividedBy(5);
        int idx = 0;
        double ae = 0.0;
        double mse = 0.0;
        for (Point expected : tsDecompressed) {
            Point actual = ts.get(idx);
            if (expected.getTimestamp() != actual.getTimestamp()) continue;
            idx++;
            double diff = actual.getValue() - expected.getValue();
            ae += Math.abs(diff);
            mse += Math.abs(diff*diff);
            assertEquals(actual.getValue(), expected.getValue(), 1.1 * epsilon, "Value did not match for timestamp " + actual.getTimestamp());
        }
        assertEquals(idx, ts.size());

        return new double[]{result.getBinary().length, result.getSegmentsSize(), (long) ae, mse};
    }


    private void runAll(String[] filenames, double[] steps) throws Exception {
        for (String filename : filenames) {
            System.out.println(filename);
            String delimiter = ",";
            TimeSeries ts = TimeSeriesReader.getTimeSeries(TestPLA.class.getResourceAsStream(filename), delimiter, true);
            for (double epsilonPct : steps) {
                long dur = 0, dedur = 0;
                for (int i=0;i<10;i++){
                    double[] simpiece = Encoding(ts.data, ts.range * epsilonPct, false, false, 1, 0.0);
                    dur += compressDuration.toNanos();
                    dedur += decompressDuration.toNanos();
                    if(i==9) {
						System.out.printf("Sim-Piece\tEpsilon: %.2f%%\tCompression Ratio: %.3f\tSegments: %d\tMAE: %.10f\tMAE%%: %.10f\tRMSE: %.10f\tRMSE%%: %.10f\tAmortizedCompressionTime: %.10f\tAmortizedDecompressionTime: %.10f\n", epsilonPct * 100, (double) ts.size / simpiece[0], (long) simpiece[1], simpiece[2]/ts.data.size(), simpiece[2]/(ts.range * ts.data.size()), Math.sqrt(simpiece[3]/ts.data.size()), Math.sqrt(simpiece[3]/(ts.data.size()))/ts.range, ((double)dur/10)/ts.data.size(), ((double)dedur/10)/ts.data.size());
					}
                }
                dur = 0;
                dedur = 0;
                for (int i=0;i<10;i++){
                    double[] simpiece = MixPiece(ts.data, ts.range * epsilonPct, false, false);
                    dur += compressDuration.toNanos();
                    dedur += decompressDuration.toNanos();
                    if(i==9) System.out.printf("Mix-Piece\tEpsilon: %.2f%%\tCompression Ratio: %.3f\tSegments: %d\tMAE: %.10f\tMAE%%: %.10f\tRMSE: %.10f\tRMSE%%: %.10f\tAmortizedCompressionTime: %.10f\tAmortizedDecompressionTime: %.10f\n", epsilonPct * 100, (double) ts.size / simpiece[0], (long) simpiece[1], simpiece[2]/ts.data.size(), simpiece[2]/(ts.range * ts.data.size()), Math.sqrt(simpiece[3]/ts.data.size()), Math.sqrt(simpiece[3]/(ts.data.size()))/ts.range, ((double)dur/10)/ts.data.size(), ((double)dedur/10)/ts.data.size());
                }
                double step = 0.05;
                for (double p = 0.0; p<1.0; p+=step){
                    dur = 0;
                    dedur = 0;
                    int iter = 3;
                    for (int i=0;i<iter;i++){
                        double[] greedy = MixPieceTunablePeekAhead(ts.data, ts.range * epsilonPct, false, false, p);
                        dur += compressDuration.toNanos();
                        dedur += decompressDuration.toNanos();
                        if(i==(iter-1))System.out.printf("TailorPieceGD(^%.2f)\tEpsilon: %.2f%%\tCompression Ratio: %.3f\tSegments: %d\tMAE: %.10f\tMAE%%: %.10f\tRMSE: %.10f\tRMSE%%: %.10f\tAmortizedCompressionTime: %.10f\tAmortizedDecompressionTime: %.10f\n", p, epsilonPct * 100, (double) ts.size / greedy[0], (long)greedy[1], greedy[2]/ts.data.size(), greedy[2]/(ts.range * ts.data.size()), Math.sqrt(greedy[3]/ts.data.size()), Math.sqrt(greedy[3]/(ts.data.size()))/ts.range, ((double)dur/iter)/ts.data.size(), ((double)dedur/iter)/ts.data.size());
                    }
                    if(p > 0.9) step = 0.01;
                }
                double[] best0 = MixPieceQuantOptimal(ts.data, ts.range * epsilonPct, false, false, 0.0);
                System.out.printf("Min-Segments\tEpsilon: %.2f%%\tCompression Ratio: %.3f\tSegments: %d\tMAE: %.10f\tMAE%%: %.10f\tRMSE: %.10f\tRMSE%%: %.10f\tAmortizedCompressionTime: %.10f\tAmortizedDecompressionTime: %.10f\n", epsilonPct * 100, (double) ts.size / best0[0], (long)best0[1], best0[2]/ts.data.size(), best0[2]/(ts.range * ts.data.size()), Math.sqrt(best0[3]/ts.data.size()), Math.sqrt(best0[3]/(ts.data.size()))/ts.range, (double) compressDuration.toNanos()/ts.data.size(), (double) decompressDuration.toNanos()/ts.data.size());
                int[] powers = {0, -1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -15, -20, -25, -30, -35, -40, -45, -50, -55, -60, -70, -80};
                for (int pow : powers) {
                    double[] bestDP = MixPieceQuantOptimal(ts.data, ts.range * epsilonPct, false, false, Math.pow(2, pow));
                    System.out.printf("TailorPieceDP(^%.14f)\tEpsilon: %.2f%%\tCompression Ratio: %.3f\tSegments: %d\tMAE: %.10f\tMAE%%: %.10f\tRMSE: %.10f\tRMSE%%: %.10f\tAmortizedCompressionTime: %.10f\tAmortizedDecompressionTime: %.10f\n", Math.pow(2, pow), epsilonPct * 100, (double) ts.size / bestDP[0], (long)bestDP[1], bestDP[2]/ts.data.size(), bestDP[2]/(ts.range * ts.data.size()), Math.sqrt(bestDP[3]/ts.data.size()), Math.sqrt(bestDP[3]/(ts.data.size()))/ts.range, (double) compressDuration.toNanos()/ts.data.size(), (double) decompressDuration.toNanos()/ts.data.size());
                }
                double[] bestDP = MixPieceMinimumAngle(ts.data, ts.range * epsilonPct, false, false);
                System.out.printf("TailorPieceMA\tEpsilon: %.2f%%\tCompression Ratio: %.3f\tSegments: %d\tMAE: %.10f\tMAE%%: %.10f\tRMSE: %.10f\tRMSE%%: %.10f\tAmortizedCompressionTime: %.10f\tAmortizedDecompressionTime: %.10f\n", epsilonPct * 100, (double) ts.size / bestDP[0], (long)bestDP[1], bestDP[2]/ts.data.size(), bestDP[2]/(ts.range * ts.data.size()), Math.sqrt(bestDP[3]/ts.data.size()), Math.sqrt(bestDP[3]/(ts.data.size()))/ts.range, (double) compressDuration.toNanos()/ts.data.size(), (double) decompressDuration.toNanos()/ts.data.size());
                System.out.println();
            }
            System.out.println();

        }
    }

    @Test
    public void TestAll() throws Exception {
         double[] steps = {0.01, 0.05, 0.1};

        String[] filenames = { "/Yoga.csv.gz", "/StarLightCurves-sample.csv.gz", "/Car.csv.gz" };

        runAll(filenames, steps);
    }
}
