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

        long compressedSize = binary.length;

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

    private double[] MixPiece2DDP(List<Point> ts, double epsilon, boolean variableByte, boolean zstd, double pow) throws Exception {
        compressDuration = Duration.ZERO;
        Instant start = Instant.now();
        Result result = MixPiece.compress2dDp(ts, epsilon, variableByte, zstd, pow);
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



    private void run(String[] filenames, double epsilonStart, double epsilonStep, double epsilonEnd) throws Exception {
        for (String filename : filenames) {
            System.out.println(filename);
            String delimiter = ",";
            TimeSeries ts = TimeSeriesReader.getTimeSeries(getClass().getResourceAsStream(filename), delimiter, true);
            for (double epsilonPct = epsilonStart; epsilonPct <= epsilonEnd; epsilonPct += epsilonStep) {
        		long dur = 0;
                /*for (int i=0;i<10;i++){
                            double[] simpiece = Encoding(ts.data, ts.range * epsilonPct, false, false, 1, 0.0);
                            dur += compressDuration.toMillis();
                            if(i==9) System.out.printf("Sim-Piece\tEpsilon: %.2f%%\tCompression Ratio: %.3f\tExecution Time: %dms\tSegments: %d\tMAE: %.10f\tMAE%%: %.10f\t%.2f\n", epsilonPct * 100, (double) ts.size / simpiece[0], dur/10, (long) simpiece[1], simpiece[2]/ts.data.size(), simpiece[2]/(ts.range * ts.data.size()), ts.range);
                }*/
                dur = 0;
                for (int i=0;i<10;i++){
                            double[] simpiece = Encoding(ts.data, ts.range * epsilonPct, false, false, 5, 0.0);
                            dur += compressDuration.toMillis();
                            if(i==9) System.out.printf("Sim-Piece*\tEpsilon: %.2f%%\tCompression Ratio: %.3f\tExecution Time: %dms\tSegments: %d\tMAE: %.10f\tMAE%%: %.10f\t%.2f\n", epsilonPct * 100, (double) ts.size / simpiece[0], dur/10, (long) simpiece[1], simpiece[2]/ts.data.size(), simpiece[2]/(ts.range * ts.data.size()), ts.range);
                }
                dur = 0;
                for (int i=0;i<10;i++){
                    double[] simpiece = MixPiece(ts.data, ts.range * epsilonPct, false, false);
                    dur += compressDuration.toMillis();
                    if(i==9) System.out.printf("Mix-Piece\tEpsilon: %.2f%%\tCompression Ratio: %.3f\tExecution Time: %dms\tSegments: %d\tMAE: %.10f\tMAE%%: %.10f\t%.2f\n", epsilonPct * 100, (double) ts.size / simpiece[0], dur/10, (long) simpiece[1], simpiece[2]/ts.data.size(), simpiece[2]/(ts.range * ts.data.size()), ts.range);
                }
                double step = 0.05;
                for (double p = 0.0; p<1.0; p+=step){
			        dur =0;
                	for (int i=0;i<5;i++){
        	        	double[] greedy = Encoding(ts.data, ts.range * epsilonPct, false, false, 2, p);
				        dur += compressDuration.toMillis();
	                	if(i==4)System.out.printf("TailorPieceGD(^%.2f)\tEpsilon: %.2f%%\tCompression Ratio: %.3f\tExecution Time: %dms\tSegments: %d\tMAE: %.10f\tMAE%%: %.10f\n", p, epsilonPct * 100, (double) ts.size / greedy[0], dur/5, (long)greedy[1], greedy[2]/ts.data.size(), greedy[2]/(ts.range * ts.data.size()));
			        }
                    if(p > 0.9) step = 0.01;
		        }
                double[] best0 = Encoding(ts.data, ts.range * epsilonPct, false, false, 0, 0.0);
                System.out.printf("Min-Segments\tEpsilon: %.2f%%\tCompression Ratio: %.3f\tExecution Time: %dms\tSegments: %d\tMAE: %.10f\tMAE%%: %.10f\n", epsilonPct * 100, (double) ts.size / best0[0], compressDuration.toMillis(), (long)best0[1], best0[2]/ts.data.size(), best0[2]/(ts.range * ts.data.size()));
                best0 = Encoding(ts.data, ts.range * epsilonPct, false, false, 0, Math.pow(2, -20));
                System.out.printf("TailorPieceDP(^%.14f)\tEpsilon: %.2f%%\tCompression Ratio: %.3f\tExecution Time: %dms\tSegments: %d\tMAE: %.10f\tMAE%%: %.10f\n", Math.pow(2, -18), epsilonPct * 100, (double) ts.size / best0[0], compressDuration.toMillis(), (long)best0[1], best0[2]/ts.data.size(), best0[2]/(ts.range * ts.data.size()));

                System.out.println();
            }
            System.out.println();

        }
    }


    private void runAll(String[] filenames, double[] steps) throws Exception {
        for (String filename : filenames) {
            System.out.println(filename);
            String delimiter = ",";
            TimeSeries ts = TimeSeriesReader.getTimeSeries(new FileInputStream(new File(filename)), delimiter, true);
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
//                double[] best2DDP = MixPiece2DDP(ts.data, ts.range * epsilonPct, false, false, Math.pow(2, pow));
//                System.out.printf("TailorPiece2D-DP(^%.8f)\tEpsilon: %.2f%%\tCompression Ratio: %.3f\tSegments: %d\tMAE: %.10f\tMAE%%: %.10f\tRMSE: %.10f\tRMSE%%: %.10f\tAmortizedCompressionTime: %.10f\tAmortizedDecompressionTime: %.10f\n", Math.pow(2, pow), epsilonPct * 100, (double) ts.size / best2DDP[0], (long)best2DDP[1], best2DDP[2]/ts.data.size(), best2DDP[2]/(ts.range * ts.data.size()), Math.sqrt(best2DDP[3]/ts.data.size()), Math.sqrt(best2DDP[3]/(ts.data.size()))/ts.range, (double) compressDuration.toNanos()/ts.data.size(), (double) decompressDuration.toNanos()/ts.data.size());

                System.out.println();
            }
            System.out.println();

        }
    }

    @Test
    public void TestAll() throws Exception {
        double[] steps = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1};
        // double[] steps = {0.01, 0.05, 0.1};
        String[] filenames = {
/*//"/home/panagiotis/workspace/cameo/data/pla-humidity.csv.gz",
//"/home/panagiotis/workspace/cameo/data/pla-irbiotemp.csv.gz",*/
//"/home/panagiotis/workspace/datasets/UCR/processed-1M/ACSF1.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/Adiac.csv.gz",
//"/home/panagiotis/workspace/datasets/UCR/processed-1M/AllGestureWiimoteX.csv.gz",
//"/home/panagiotis/workspace/datasets/UCR/processed-1M/AllGestureWiimoteY.csv.gz",
//"/home/panagiotis/workspace/datasets/UCR/processed-1M/AllGestureWiimoteZ.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/ArrowHead.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/Beef.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/BeetleFly.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/BirdChicken.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/BME.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/Car.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/CBF.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/Chinatown.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/ChlorineConcentration.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/CinCECGTorso-sample.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/Coffee.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/Computers-sample.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/CricketX.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/CricketY.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/CricketZ.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/Crop.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/DiatomSizeReduction.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/DistalPhalanxOutlineAgeGroup.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/DistalPhalanxOutlineCorrect.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/DistalPhalanxTW.csv.gz",
//"/home/panagiotis/workspace/datasets/UCR/processed-1M/DodgerLoopDay.csv.gz",
//"/home/panagiotis/workspace/datasets/UCR/processed-1M/DodgerLoopGame.csv.gz",
//"/home/panagiotis/workspace/datasets/UCR/processed-1M/DodgerLoopWeekend.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/Earthquakes.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/ECG200.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/ECG5000.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/ECGFiveDays.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/ElectricDevices.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/EOGHorizontalSignal-sample.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/EOGVerticalSignal-sample.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/EthanolLevel-sample.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/FaceAll.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/FaceFour.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/FacesUCR.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/FiftyWords.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/Fish.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/FordA.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/FordB.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/FreezerRegularTrain.csv.gz",
//"/home/panagiotis/workspace/datasets/UCR/processed-1M/FreezerSmallTrain.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/Fungi.csv.gz",
//"/home/panagiotis/workspace/datasets/UCR/processed-1M/GestureMidAirD1.csv.gz",
//"/home/panagiotis/workspace/datasets/UCR/processed-1M/GestureMidAirD2.csv.gz",
//"/home/panagiotis/workspace/datasets/UCR/processed-1M/GestureMidAirD3.csv.gz",
//"/home/panagiotis/workspace/datasets/UCR/processed-1M/GesturePebbleZ1.csv.gz",
//"/home/panagiotis/workspace/datasets/UCR/processed-1M/GesturePebbleZ2.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/GunPointAgeSpan.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/GunPoint.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/GunPointMaleVersusFemale.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/GunPointOldVersusYoung.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/Ham.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/HandOutlines-sample.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/Haptics-sample.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/Herring.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/HouseTwenty-sample.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/InlineSkate-sample.csv.gz",
//"/home/panagiotis/workspace/datasets/UCR/processed-1M/InsectEPGRegularTrain.csv.gz",
//"/home/panagiotis/workspace/datasets/UCR/processed-1M/InsectEPGSmallTrain.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/InsectWingbeatSound.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/ItalyPowerDemand.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/LargeKitchenAppliances-sample.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/Lightning2.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/Lightning7.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/Mallat.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/Meat.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/MedicalImages.csv.gz",
//"/home/panagiotis/workspace/datasets/UCR/processed-1M/MelbournePedestrian-sample.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/MiddlePhalanxOutlineAgeGroup.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/MiddlePhalanxOutlineCorrect.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/MiddlePhalanxTW.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/MixedShapesRegularTrain.csv.gz",
//"/home/panagiotis/workspace/datasets/UCR/processed-1M/MixedShapesSmallTrain.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/MoteStrain.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/NonInvasiveFetalECGThorax1.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/NonInvasiveFetalECGThorax2.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/OliveOil.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/OSULeaf.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/PhalangesOutlinesCorrect.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/Phoneme.csv.gz",
//"/home/panagiotis/workspace/datasets/UCR/processed-1M/PickupGestureWiimoteZ.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/PigAirwayPressure-sample.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/PigArtPressure.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/PigCVP-sample.csv.gz",
//"/home/panagiotis/workspace/datasets/UCR/processed-1M/PLAID.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/Plane.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/PowerCons.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/ProximalPhalanxOutlineAgeGroup.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/ProximalPhalanxOutlineCorrect.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/ProximalPhalanxTW.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/RefrigerationDevices.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/Rock-sample.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/ScreenType-sample.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/SemgHandGenderCh2.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/SemgHandMovementCh2.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/SemgHandSubjectCh2.csv.gz",
//"/home/panagiotis/workspace/datasets/UCR/processed-1M/ShakeGestureWiimoteZ.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/ShapeletSim.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/ShapesAll.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/SmallKitchenAppliances-sample.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/SmoothSubspace.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/SonyAIBORobotSurface1.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/SonyAIBORobotSurface2.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/StarLightCurves.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/Strawberry.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/SwedishLeaf.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/Symbols.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/SyntheticControl.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/ToeSegmentation1.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/ToeSegmentation2.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/Trace.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/TwoLeadECG.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/TwoPatterns.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/UMD.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/UWaveGestureLibraryAll.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/UWaveGestureLibraryX.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/UWaveGestureLibraryY.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/UWaveGestureLibraryZ.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/Wafer.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/Wine.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/WordSynonyms.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/Worms.csv.gz",
//"/home/panagiotis/workspace/datasets/UCR/processed-1M/WormsTwoClass.csv.gz",
"/home/panagiotis/workspace/datasets/UCR/processed-1M/Yoga.csv.gz"
	};

        //String[] filenames = {"/Yoga.csv.gz", "/Rock.csv.gz", "/Worms.csv.gz", "/Trace.csv.gz", "/StarLightCurves-sample.csv.gz", "/Car.csv.gz", "/CinCECGTorso-sample2.csv.gz", "/Plane.csv.gz", "/citytemp_f32_sample.csv.gz", "/jane_street_f64_sample.csv.gz", "/solar_wind_f32_sample.csv.gz", "/Lightning.csv.gz", "/Cricket.csv.gz", "/FaceFour.csv.gz", "/WindSpeed_sample.csv.gz" };
//        String[] filenames = { "/Yoga-sample.csv.gz", "/Yoga-sample-2.csv.gz", "/Yoga-sample-3.csv.gz", "/Yoga-sample-4.csv.gz", "/Yoga-sample-5.csv.gz", "/StarLightCurves-sample.csv.gz", "/Car.csv.gz" };
//        String[] filenames = { "/Yoga.csv.gz", "/StarLightCurves-sample.csv.gz", "/Car.csv.gz" };
//        String[] filenames = {"/Yoga.csv.gz" };

        runAll(filenames, steps);
    }
}
