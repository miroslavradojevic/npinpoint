package auxiliary;

import ij.ImagePlus;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

/**
 * Automatic thresholding technique based on the entopy of the histogram.
 * See: P.K. Sahoo, S. Soltani, K.C. Wong and, Y.C. Chen "A Survey of Thresholding Techniques", Computer Vision, Graphics, and Image Processing, Vol. 41, pp.233-260, 1988.
 */
public class EntrTh implements PlugInFilter {

    public int setup(String s, ImagePlus imagePlus) {
        return PlugInFilter.DOES_8G | PlugInFilter.DOES_STACKS;
    }

    public void run(ImageProcessor imageProcessor) {
        int[] hist = imageProcessor.getHistogram();
        int threshold = entropySplit(hist);
        System.out.println("thr: " + threshold);
        imageProcessor.threshold(threshold);
    }

    public ByteProcessor runMaxEntropyThreshold(ImageProcessor _in_float_processor) {

        ByteProcessor bp = _in_float_processor.convertToByteProcessor();
        int[] hist = bp.getHistogram();
        int threshold = entropySplit(hist);
        bp.threshold(threshold);

        return bp;

    }

    public ByteProcessor runMaxEntropyThreshold(ImageProcessor _in_float_processor, float _detection_sensitivity) {

        ByteProcessor bp = _in_float_processor.convertToByteProcessor();
        int[] hist = bp.getHistogram();
        int thresholdMaxEntropy = entropySplit(hist);

        int thresholdMax1 = Integer.MIN_VALUE;
        int thresholdMin1 = Integer.MAX_VALUE;
        byte[] vals = (byte[]) bp.getPixels();
        for (int i = 0; i < vals.length; i++) {

            int v = vals[i]&0xff;

            if (v>thresholdMax1) thresholdMax1 = v;
            if (v<thresholdMin1) thresholdMin1 = v;

        }

        // choose threshold depending on the sensitivity
        int threshold = (int) (_detection_sensitivity * (thresholdMaxEntropy/2f) + (1 - _detection_sensitivity) * thresholdMax1); // thresholdMaxEntropy thresholdMin1

        bp.threshold(threshold);

        return bp;

    }

    private int entropySplit(int[] hist) {

        // Normalize histogram, that is makes the sum of all bins equal to 1.
        double sum = 0;
        for (int i = 0; i < hist.length; ++i) {
            sum += hist[i];
        }
        if (sum == 0) {
            // This should not normally happen, but...
            throw new IllegalArgumentException("Empty histogram: sum of all bins is zero.");
        }

        double[] normalizedHist = new double[hist.length];
        for (int i = 0; i < hist.length; i++) {
            normalizedHist[i] = hist[i] / sum;
        }

        //
        double[] pT = new double[hist.length];
        pT[0] = normalizedHist[0];
        for (int i = 1; i < hist.length; i++) {
            pT[i] = pT[i - 1] + normalizedHist[i];
        }

        // Entropy for black and white parts of the histogram
        final double epsilon = Double.MIN_VALUE;
        double[] hB = new double[hist.length];
        double[] hW = new double[hist.length];
        for (int t = 0; t < hist.length; t++) {
            // Black entropy
            if (pT[t] > epsilon) {
                double hhB = 0;
                for (int i = 0; i <= t; i++) {
                    if (normalizedHist[i] > epsilon) {
                        hhB -= normalizedHist[i] / pT[t] * Math.log(normalizedHist[i] / pT[t]);
                    }
                }
                hB[t] = hhB;
            } else {
                hB[t] = 0;
            }

            // White  entropy
            double pTW = 1 - pT[t];
            if (pTW > epsilon) {
                double hhW = 0;
                for (int i = t + 1; i < hist.length; ++i) {
                    if (normalizedHist[i] > epsilon) {
                        hhW -= normalizedHist[i] / pTW * Math.log(normalizedHist[i] / pTW);
                    }
                }
                hW[t] = hhW;
            } else {
                hW[t] = 0;
            }
        }

        // Find histogram index with maximum entropy
        double jMax = hB[0] + hW[0];
        int tMax = 0;
        for (int t = 1; t < hist.length; ++t) {
            double j = hB[t] + hW[t];
            if (j > jMax) {
                jMax = j;
                tMax = t;
            }
        }

        return tMax;
    }

}
