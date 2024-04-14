package org.pla.compression.benchmarks;

import org.pla.compression.pmcmr.PMCMR;
import org.pla.compression.simpiece.SimPiece;
import org.pla.compression.swingfilter.SwingFilter;
import org.pla.compression.util.Point;
import org.pla.compression.util.TimeSeries;
import org.pla.compression.util.TimeSeriesReader;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.time.Duration;
import java.time.Instant;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class TestPLA {
    private Duration duration;

    private long PMCMR(List<Point> ts, double epsilon) throws IOException {
        duration = Duration.ZERO;
        Instant start = Instant.now();
        PMCMR pmcmr = new PMCMR(ts, epsilon);
        duration = Duration.between(start, Instant.now());

        byte[] binary = pmcmr.toByteArray();

        long compressedSize = binary.length;

        pmcmr = new PMCMR(binary);
        List<Point> tsDecompressed = pmcmr.decompress();
        int idx = 0;
        for (Point expected : tsDecompressed) {
            Point actual = ts.get(idx);
            if (expected.getTimestamp() != actual.getTimestamp()) continue;
            idx++;
            assertEquals(actual.getValue(), expected.getValue(), 1.1 * epsilon, "Value did not match for timestamp " + actual.getTimestamp());
        }
        assertEquals(idx, ts.size());

        return compressedSize;
    }


    private long Swing(List<Point> ts, double epsilon) throws IOException {
        duration = Duration.ZERO;
        Instant start = Instant.now();
        SwingFilter swingFilter = new SwingFilter(ts, epsilon);
        duration = Duration.between(start, Instant.now());
        byte[] binary = swingFilter.toByteArray();

        long compressedSize = binary.length;

        swingFilter = new SwingFilter(binary);
        List<Point> tsDecompressed = swingFilter.decompress();
        int idx = 0;
        for (Point expected : tsDecompressed) {
            Point actual = ts.get(idx);
            if (expected.getTimestamp() != actual.getTimestamp()) continue;
            idx++;
            assertEquals(actual.getValue(), expected.getValue(), 1.1 * epsilon, "Value did not match for timestamp " + actual.getTimestamp());
        }
        assertEquals(idx, ts.size());

        return compressedSize;
    }


    private double[] SimPiece(List<Point> ts, double epsilon, boolean variableByte, boolean zstd, int mode, double pow) throws IOException {
        duration = Duration.ZERO;
        Instant start = Instant.now();
        SimPiece simPiece = new SimPiece(ts, epsilon, mode, pow);
        duration = Duration.between(start, Instant.now());
        byte[] binary = simPiece.toByteArray(variableByte, zstd);

        long compressedSize = binary.length;

        simPiece = new SimPiece(binary, variableByte, zstd);
        List<Point> tsDecompressed = simPiece.decompress();
        int idx = 0;
        double ae = 0.0;
        for (Point expected : tsDecompressed) {
            Point actual = ts.get(idx);
            if (expected.getTimestamp() != actual.getTimestamp()) continue;
            idx++;
            ae += Math.abs(actual.getValue() - expected.getValue());
            assertEquals(actual.getValue(), expected.getValue(), 1.1 * epsilon, "Value did not match for timestamp " + actual.getTimestamp());
        }
        assertEquals(idx, ts.size());

        return new double[]{compressedSize, simPiece.getSegments().size(), (long) ae};
    }


    private void run(String[] filenames, double epsilonStart, double epsilonStep, double epsilonEnd) throws IOException {
        for (String filename : filenames) {
            System.out.println(filename);
            String delimiter = ",";
            TimeSeries ts = TimeSeriesReader.getTimeSeries(getClass().getResourceAsStream(filename), delimiter, true);
            for (double epsilonPct = epsilonStart; epsilonPct <= epsilonEnd; epsilonPct += epsilonStep) {
		long dur = 0;
		for (int i=0;i<10;i++){
	                double[] simpiece = SimPiece(ts.data, ts.range * epsilonPct, false, false, 1, 0.0);
            		dur += duration.toMillis();
                	if(i==9) System.out.printf("Sim-Piece\tEpsilon: %.2f%%\tCompression Ratio: %.3f\tExecution Time: %dms\tSegments: %d\tMAE: %.10f\tMAE%%: %.10f\t%.2f\n", epsilonPct * 100, (double) ts.size / simpiece[0], dur/10, (long) simpiece[1], simpiece[2]/ts.data.size(), simpiece[2]/(ts.range * ts.data.size()), ts.range);
		}
//                for (double pow=0.0005; pow<=5; pow*=2) {
//                    long[] greedy2 = SimPiece(ts.data, ts.range * epsilonPct, false, false, 4, pow);
//                    System.out.printf("2nd-Variation (^%.8f)\tEpsilon: %.2f%%\tCompression Ratio: %.3f\tExecution Time: %dms\tSegments: %d\n", pow, epsilonPct * 100, (double) ts.size / greedy2[0], duration.toMillis(), greedy2[1]);
//                }
                for (double p = 0.0; p<1.00; p+=0.05){
			dur =0;
                	for (int i=0;i<5;i++){
        	        	double[] greedy = SimPiece(ts.data, ts.range * epsilonPct, false, false, 2, p);
				dur += duration.toMillis();	
	                	if(i==4)System.out.printf("Two-step Greedy(^%.2f)\tEpsilon: %.2f%%\tCompression Ratio: %.3f\tExecution Time: %dms\tSegments: %d\tMAE: %.10f\tMAE%%: %.10f\n", p, epsilonPct * 100, (double) ts.size / greedy[0], dur/5, (long)greedy[1], greedy[2]/ts.data.size(), greedy[2]/(ts.range * ts.data.size()));
			}
                    break;
		}
                //greedy = SimPiece(ts.data, ts.range * epsilonPct, false, false, 2, 0.5);
                //System.out.printf("Two-step Greedy/0.5\tEpsilon: %.2f%%\tCompression Ratio: %.3f\tExecution Time: %dms\tSegments: %d\n", epsilonPct * 100, (double) ts.size / greedy[0], duration.toMillis(), greedy[1]);
                //greedy = SimPiece(ts.data, ts.range * epsilonPct, false, false, 2, 0.7);
                //System.out.printf("Two-step Greedy/0.7\tEpsilon: %.2f%%\tCompression Ratio: %.3f\tExecution Time: %dms\tSegments: %d\n", epsilonPct * 100, (double) ts.size / greedy[0], duration.toMillis(), greedy[1]);
                //greedy = SimPiece(ts.data, ts.range * epsilonPct, false, false, 2, 0.9);
                //System.out.printf("Two-step Greedy/0.9\tEpsilon: %.2f%%\tCompression Ratio: %.3f\tExecution Time: %dms\tSegments: %d\n", epsilonPct * 100, (double) ts.size / greedy[0], duration.toMillis(), greedy[1]);
                /*long[] greedy2 = SimPiece(ts.data, ts.range * epsilonPct, false, false, 3, 2);
                System.out.printf("Two-step Greedy^\tEpsilon: %.2f%%\tCompression Ratio: %.3f\tExecution Time: %dms\tSegments: %d\n", epsilonPct * 100, (double) ts.size / greedy2[0], duration.toMillis(), greedy2[1]);*/

                for (double pow=0.0, i=-18; pow<=2; pow=Math.pow(2, i++)) {
                    double[] best0 = SimPiece(ts.data, ts.range * epsilonPct, false, false, 0, pow);
                    System.out.printf("Optimal (^%.8f)\tEpsilon: %.2f%%\tCompression Ratio: %.3f\tExecution Time: %dms\tSegments: %d\tMAE: %.10f\tMAE%%: %.10f\n", pow, epsilonPct * 100, (double) ts.size / best0[0], duration.toMillis(), (long)best0[1], best0[2]/ts.data.size(), best0[2]/(ts.range * ts.data.size()));
			    //if(pow>0.0)
                break;
                }
                System.out.println();
            }
            System.out.println();

//            System.out.println("Sim-Piece Variable Byte");
//            for (double epsilonPct = epsilonStart; epsilonPct <= epsilonEnd; epsilonPct += epsilonStep)
//                System.out.printf("Epsilon: %.2f%%\tCompression Ratio: %.3f\tExecution Time: %dms\n", epsilonPct * 100, (double) ts.size / SimPiece(ts.data, ts.range * epsilonPct, true, false), duration.toMillis());
//
//
//            System.out.println("Sim-Piece Variable Byte & ZStd");
//            for (double epsilonPct = epsilonStart; epsilonPct <= epsilonEnd; epsilonPct += epsilonStep)
//                System.out.printf("Epsilon: %.2f%%\tCompression Ratio: %.3f\tExecution Time: %dms\n", epsilonPct * 100, (double) ts.size / SimPiece(ts.data, ts.range * epsilonPct, true, true), duration.toMillis());
//
//            System.out.println("Swing");
//            for (double epsilonPct = epsilonStart; epsilonPct <= epsilonEnd; epsilonPct += epsilonStep)
//                System.out.printf("Epsilon: %.2f%%\tCompression Ratio: %.3f\tExecution Time: %dms\n", epsilonPct * 100, (double) ts.size / Swing(ts.data, ts.range * epsilonPct), duration.toMillis());
//
//            System.out.println("PMCMR");
//            for (double epsilonPct = epsilonStart; epsilonPct <= epsilonEnd; epsilonPct += epsilonStep)
//                System.out.printf("Epsilon: %.2f%%\tCompression Ratio: %.3f\tExecution Time: %dms\n", epsilonPct * 100, (double) ts.size / PMCMR(ts.data, ts.range * epsilonPct), duration.toMillis());

            System.out.println();
        }
    }


    @Test
    public void TestCRAndTime() throws IOException {
        double epsilonStart = 0.005;
        double epsilonStep = 0.0025;
        double epsilonEnd = 0.051;

        String[] filenames = {"/citytemp_f32_sample.csv.gz", "/jane_street_f64_sample.csv.gz", "/nyc_taxi2015_f64_sample.csv.gz", /*"/phone_gyro_f64_sample.csv.gz",*/ "/solar_wind_f32_sample.csv.gz", /*"/ts_gas_f32_sample.csv.gz", "/wesad_chest_f64_sample.csv.gz",*/ "/Lightning.csv.gz", "/spain_gas_price_f64_sample.csv.gz", "/Cricket.csv.gz", "/FaceFour.csv.gz", "/MoteStrain.csv.gz", "/Wafer.csv.gz", "/WindSpeed_sample.csv.gz", "/WindDirection.csv.gz", "/Pressure_sample.csv.gz" };
        run(filenames, epsilonStart, epsilonStep, epsilonEnd);
    }
}
