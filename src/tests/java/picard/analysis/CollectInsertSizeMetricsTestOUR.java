package picard.analysis;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import static org.testng.Assert.fail;
import static picard.analysis.SinglePassSamProgram.POISON_PILL;
import static picard.analysis.SinglePassSamProgram.POISON_PILL_TAG;

/**
 * Created by student on 6/30/16.
 */
public class CollectInsertSizeMetricsTestOUR  extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR_SAM = new File("testdata/picard/sam");
    private static final File TEST_DATA_DIR_OUTPUT = new File("testdata/picard/analysis/our_output");
    private static final File TEST_DATA_DIR_BAM = new File("testdata/picard/analysis/OurTestFiles/");


    @Override
    public String getCommandLineProgramName() {
        return CollectInsertSizeMetrics.class.getSimpleName();
    }

    @Test
    public void testSam() throws IOException {


        final File input = new File(TEST_DATA_DIR_SAM, "insert_size_metrics_test.sam");
        final File outfile   = new File(TEST_DATA_DIR_OUTPUT, "insert_size_metrics_test.sam.insert_size_metrics");
        final File pdf   = new File(TEST_DATA_DIR_OUTPUT, "insert_size_metrics_test.sam.pdf");

        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "HISTOGRAM_FILE=" + pdf.getAbsolutePath()
        };


        //start measure worktime
        long startTime = System.nanoTime();

        runPicardCommandLine(args);

        //end measure worktime ant print it
        long endTime = System.nanoTime();
        double estTime = ((endTime-startTime)/(Math.pow(10, 9)));
        double finalValue = Math.round( estTime * 1000.0 ) / 1000.0;
        System.out.print(finalValue + "\t");
    }

    @Test
    public void testHistogramWidthIsSetProperly() throws IOException {
        final File input = new File(TEST_DATA_DIR_SAM, "insert_size_metrics_test.sam");
        final File outfile = new File(TEST_DATA_DIR_OUTPUT, "test_hist.insert_size_metrics");
        final File pdf = new File(TEST_DATA_DIR_OUTPUT, "test_hist.pdf");

        final String[] args = new String[] {
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "HISTOGRAM_FILE=" + pdf.getAbsolutePath(),
                "LEVEL=null",
                "LEVEL=READ_GROUP"
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);
        final MetricsFile<InsertSizeMetrics, Comparable<?>> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        Assert.assertEquals(output.getAllHistograms().size(), 5);
    }


    @Test
    public void testBamFromTheNet() throws IOException {

        final File input = new File(TEST_DATA_DIR_BAM, "mapt.NA12156.altex.bam");
        final File outfile   = new File(TEST_DATA_DIR_OUTPUT, "mapt.NA12156.altex.bam.insert_size_metrics");
        final File pdf   = new File(TEST_DATA_DIR_OUTPUT, "mapt.NA12156.altex.bam.pdf");

        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "HISTOGRAM_FILE=" + pdf.getAbsolutePath(),
                "VALIDATION_STRINGENCY=LENIENT"

        };


        //start measure worktime
        long startTime = System.nanoTime();

        runPicardCommandLine(args);

        //end measure worktime ant print it
        long endTime = System.nanoTime();
        double estTime = ((endTime-startTime)/(Math.pow(10, 9)));
        double finalValue = Math.round( estTime * 1000.0 ) / 1000.0;
        System.out.print(finalValue + "\t");
    }

    @Test
    public void testPoisonPill() throws Exception {
        SAMRecord record = POISON_PILL;

        if (!record.getReadPairedFlag() ||
                record.getReadUnmappedFlag() ||
                record.getMateUnmappedFlag() ||
                record.getFirstOfPairFlag() ||
                record.isSecondaryOrSupplementary() ||
                (record.getDuplicateReadFlag()) ||
                record.getInferredInsertSize() == 0 && record.getAttribute(POISON_PILL_TAG) == null) {
            fail();
        }
        System.out.println(record.getAttribute(POISON_PILL_TAG));
    }

    @Test
    public void testBigBamFromTheNet() throws IOException {

        final File input = new File(TEST_DATA_DIR_BAM, "HG00117.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam");
        final File outfile   = new File(TEST_DATA_DIR_OUTPUT, "HG00117.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam.insert_size_metrics");
        final File pdf   = new File(TEST_DATA_DIR_OUTPUT, "HG00117.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam.pdf");

        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "HISTOGRAM_FILE=" + pdf.getAbsolutePath(),
                "VALIDATION_STRINGENCY=LENIENT"
        };


        //start measure worktime
        long startTime = System.nanoTime();

        runPicardCommandLine(args);

        //end measure worktime ant print it
        long endTime = System.nanoTime();
        double estTime = ((endTime-startTime)/(Math.pow(10, 9)));
        double finalValue = Math.round( estTime * 1000.0 ) / 1000.0;
        System.out.print(finalValue + "\t");
    }
}
