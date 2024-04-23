package org.pla.compression.benchmarks;

import org.pla.compression.encodings.Encoding;
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

    private double[] Encoding(List<Point> ts, double epsilon, boolean variableByte, boolean zstd, int mode, double pow) throws IOException {
        duration = Duration.ZERO;
        Instant start = Instant.now();
        Encoding encoding = new Encoding(ts, epsilon, mode, pow);
        duration = Duration.between(start, Instant.now());
        byte[] binary = encoding.toByteArray(variableByte, zstd);

        long compressedSize = binary.length;

        encoding = new Encoding(binary, variableByte, zstd);
        List<Point> tsDecompressed = encoding.decompress();
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

        return new double[]{compressedSize, encoding.getSegments().size(), (long) ae};
    }


    private void run(String[] filenames, double epsilonStart, double epsilonStep, double epsilonEnd) throws IOException {
        for (String filename : filenames) {
            System.out.println(filename);
            String delimiter = ",";
            TimeSeries ts = TimeSeriesReader.getTimeSeries(getClass().getResourceAsStream(filename), delimiter, true);
            for (double epsilonPct = epsilonStart; epsilonPct <= epsilonEnd; epsilonPct += epsilonStep) {
        		long dur = 0;
                for (int i=0;i<10;i++){
                            double[] simpiece = Encoding(ts.data, ts.range * epsilonPct, false, false, 1, 0.0);
                            dur += duration.toMillis();
                            if(i==9) System.out.printf("Sim-Piece\tEpsilon: %.2f%%\tCompression Ratio: %.3f\tExecution Time: %dms\tSegments: %d\tMAE: %.10f\tMAE%%: %.10f\t%.2f\n", epsilonPct * 100, (double) ts.size / simpiece[0], dur/10, (long) simpiece[1], simpiece[2]/ts.data.size(), simpiece[2]/(ts.range * ts.data.size()), ts.range);
                }
                double step = 0.05;
                for (double p = 0.0; p<1.0; p+=step){
			        dur =0;
                	for (int i=0;i<5;i++){
        	        	double[] greedy = Encoding(ts.data, ts.range * epsilonPct, false, false, 2, p);
				        dur += duration.toMillis();
	                	if(i==4)System.out.printf("Tunable-Peek-Ahead(^%.2f)\tEpsilon: %.2f%%\tCompression Ratio: %.3f\tExecution Time: %dms\tSegments: %d\tMAE: %.10f\tMAE%%: %.10f\n", p, epsilonPct * 100, (double) ts.size / greedy[0], dur/5, (long)greedy[1], greedy[2]/ts.data.size(), greedy[2]/(ts.range * ts.data.size()));
			        }
                    if(p > 0.9) step = 0.01;
		        }
                double[] best0 = Encoding(ts.data, ts.range * epsilonPct, false, false, 0, 0.0);
                System.out.printf("Quant-Optimal\tEpsilon: %.2f%%\tCompression Ratio: %.3f\tExecution Time: %dms\tSegments: %d\tMAE: %.10f\tMAE%%: %.10f\n", epsilonPct * 100, (double) ts.size / best0[0], duration.toMillis(), (long)best0[1], best0[2]/ts.data.size(), best0[2]/(ts.range * ts.data.size()));
                best0 = Encoding(ts.data, ts.range * epsilonPct, false, false, 0, Math.pow(2, -18));
                System.out.printf("Angle-Aware(^%.8f)\tEpsilon: %.2f%%\tCompression Ratio: %.3f\tExecution Time: %dms\tSegments: %d\tMAE: %.10f\tMAE%%: %.10f\n", Math.pow(2, -18), epsilonPct * 100, (double) ts.size / best0[0], duration.toMillis(), (long)best0[1], best0[2]/ts.data.size(), best0[2]/(ts.range * ts.data.size()));

                System.out.println();
            }
            System.out.println();

        }
    }


    @Test
    public void TestCRAndTime() throws IOException {
        double epsilonStart = 0.01;
        double epsilonStep = 0.005;
        double epsilonEnd = 0.051;

//        String[] filenames = {"/citytemp_f32_sample.csv.gz", "/jane_street_f64_sample.csv.gz", "/nyc_taxi2015_f64_sample.csv.gz", /*"/phone_gyro_f64_sample.csv.gz",*/ "/solar_wind_f32_sample.csv.gz", /*"/ts_gas_f32_sample.csv.gz", "/wesad_chest_f64_sample.csv.gz",*/ "/Lightning.csv.gz", "/spain_gas_price_f64_sample.csv.gz", "/Cricket.csv.gz", "/FaceFour.csv.gz", "/MoteStrain.csv.gz", "/Wafer.csv.gz", "/WindSpeed_sample.csv.gz", "/WindDirection.csv.gz" /*"/Pressure_sample.csv.gz"*/ };
	    String[] filenames = {"/jane_street_f64_sample.csv.gz"};
        run(filenames, epsilonStart, epsilonStep, epsilonEnd);
    }
}
