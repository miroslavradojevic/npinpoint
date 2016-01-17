import ij.ImageStack;
import ij.gui.Plot;
import ij.process.ImageProcessor;

/** use selected list of foreground locations and extracts profiles from those locations */
public class Profiler extends Thread {

    private int begN, endN;

    public static float[][] inimg_xy;                   // input image
    public static Sphere2D  sph2d;                      // sphere class - spherical computations
    public static int[][] 	i2xy;                       // selected locations mapping
    public static int[][]   xy2i;                       // mapping

    // output
    public static short[][]	prof2;

    public Profiler(int n0, int n1)
    {
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadTemplate(Sphere2D _sph2d, int[][] _i2xy, int[][] _xy2i, float[][] _inimg_xy){

        sph2d = _sph2d;                     // there will be one Sphere3D instance created at the beginning and
        i2xy = _i2xy;                       // list of foreground locations with corresponding x and y
        xy2i = _xy2i;                       // mapping assign
        inimg_xy = _inimg_xy;               // input image

        // allocate output
        prof2 = new short[i2xy.length][sph2d.getProfileLength()]; // output profiles

    }

    public void run() {

        for (int profileComponentIdx = begN; profileComponentIdx < endN; profileComponentIdx++) {

            for (int locIdx = 0; locIdx < i2xy.length; locIdx++) {

                int atX = i2xy[locIdx][0];
                int atY = i2xy[locIdx][1];

                prof2[locIdx][profileComponentIdx] = sph2d.extractProfile(profileComponentIdx, atX, atY, inimg_xy);

            }

        }
    }

    public static ImageStack getProfile(int atX, int atY){  // reads from prof2 array class member

        ImageStack is_out = new ImageStack(528, 255);

        int idx = xy2i[atX][atY];
        if (idx != -1) {
            int len = prof2[0].length;
            float[] f = new float[len];
            float[] fx = new float[len];

            for (int i=0; i<len; i++) {
                f[i] = ((prof2[idx][i] & 0xffff) / 65535f) * 255f; // retrieve the profile
                //fx[i] = i; // abscissa in cnt
                fx[i] = (i / (float) len) * 360; // abscissa in degs
            }

            Plot p = new Plot("profile at ("+atX+","+atY+")", "", "filtered", fx, f);
            p.setLineWidth(2);
            is_out.addSlice(p.getProcessor());
        }
        else {

            float[] fx = new float[sph2d.getProfileLength()];
            for (int i=0; i<fx.length; i++) fx[i] = ((i/(float)fx.length)*360);
            Plot p = new Plot("", "", "", fx, new float[fx.length]);
            is_out.addSlice(p.getProcessor());
        }

        return is_out;

    }

}