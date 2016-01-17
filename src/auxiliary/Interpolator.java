package auxiliary;

import ij.process.FloatProcessor;

 /** bilinear interpolation of pixel values in images or image stacks */
public final class Interpolator {

    public static final float interpolateAt(double x, double y, FloatProcessor inip2d) {

        float value = 0;

        boolean isIn =
                y>=0 && y<=(inip2d.getHeight()-1) &&
                        x>=0 && x<=(inip2d.getWidth()-1);

        if(isIn){

            int[] p11 = {(int)Math.floor(y),	(int)Math.floor(x)};
            int[] p12 = {(int)Math.floor(y), 	(int)Math.ceil(x)};
            int[] p21 = {(int)Math.ceil(y), 	(int)Math.floor(x)};
            int[] p22 = {(int)Math.ceil(y), 	(int)Math.ceil(x)};

            float I11_1 = inip2d.getf(p11[0]*inip2d.getWidth()+p11[1]);
            float I12_1 = inip2d.getf(p12[0]*inip2d.getWidth()+p12[1]);
            float I21_1 = inip2d.getf(p21[0]*inip2d.getWidth()+p21[1]);
            float I22_1 = inip2d.getf(p22[0]*inip2d.getWidth()+p22[1]);

            float a = (p12[1]!=p11[1])?((float)x-p11[1])/(p12[1]-p11[1]) : 0.5f;	//col
            float b = (p21[0]!=p11[0])?((float)y-p11[0])/(p21[0]-p11[0]) : 0.5f;	//row

            float I_1 = (1-a)*(1-b)*I11_1 + (1-a)*b*I21_1 + a*(1-b)*I12_1 + a*b*I22_1;

            value = I_1;

        }

        return value;

    }

    public static final float interpolateAt(float atX, float atY, float atZ, float[][][] img3d_xyz) {

        int x1 = (int) atX;
        int x2 = x1 + 1;
        float x_frac = atX - x1;

        int y1 = (int) atY;
        int y2 = y1 + 1;
        float y_frac = atY - y1;

        int z1 = (int) atZ;
        int z2 = z1 + 1;
        float z_frac = atZ - z1;

        boolean isIn =
                y1>=0 && y1<img3d_xyz[0].length &&
                        y2>=0 && y2<img3d_xyz[0].length &&
                        x1>=0 && x1<img3d_xyz.length &&
                        x2>=0 && x2<img3d_xyz.length &&
                        z1>=0 && z1<img3d_xyz[0][0].length &&
                        z2>=0 && z2<img3d_xyz[0][0].length;

        if(!isIn)
            System.out.println("\nthere was a problem with interpolation at " + atX + " , " + atY + "  , " + atZ + " -> x range: " + x1 + " -- " + x2 + " | y range: " + y1 + " -- " + y2 + " | z range: " + z1 + " -- " + z2 +"\n");

        // take neighbourhood
        float I11_1 = img3d_xyz[ x1 ][ y1 ][ z1 ];  // upper left
        float I12_1 = img3d_xyz[ x2 ][ y1 ][ z1 ];  // upper right
        float I21_1 = img3d_xyz[ x1 ][ y2 ][ z1 ]; // bottom left
        float I22_1 = img3d_xyz[ x2 ][ y2 ][ z1 ]; // bottom right

        float I11_2 = img3d_xyz[ x1 ][ y1 ][ z2 ]; // upper left
        float I12_2 = img3d_xyz[ x2 ][ y1 ][ z2 ]; // upper right
        float I21_2 = img3d_xyz[ x1 ][ y2 ][ z2 ]; // bottom left
        float I22_2 = img3d_xyz[ x2 ][ y2 ][ z2 ]; // bottom right

        return (1-z_frac)  * (  (1-y_frac) * ((1-x_frac)*I11_1 + x_frac*I12_1) + (y_frac) * ((1-x_frac)*I21_1 + x_frac*I22_1) )   +
                z_frac      * (  (1-y_frac) * ((1-x_frac)*I11_2 + x_frac*I12_2) + (y_frac) * ((1-x_frac)*I21_2 + x_frac*I22_2) );

    }

    public static final float interpolateAt(float atX, float atY, float[][] img2d_xy) {

        int x1 = (int) atX;
        int x2 = x1 + 1;
        float x_frac = atX - x1;

        int y1 = (int) atY;
        int y2 = y1 + 1;
        float y_frac = atY - y1;

        // technically this check should be there but is skipped to make things faster
        boolean isIn =
                y1>=0 && y1<img2d_xy[0].length &&
                        y2>=0 && y2<img2d_xy[0].length &&
                        x1>=0 && x1<img2d_xy.length &&
                        x2>=0 && x2<img2d_xy.length;

        if(!isIn)
            System.out.println("\nthere was a problem with interpolation at " + atX + " , " + atY + "   " + " -> x range: " + x1 + " -- " + x2 + " | y range: " + y1 + " -- " + y2 + ", " + "  \n");

        // take neighbourhood
        float I11_1 = img2d_xy[ x1  ][ y1 ];  // upper left
        float I12_1 = img2d_xy[ x2  ][ y1 ];  // upper right
        float I21_1 = img2d_xy[ x1  ][ y2  ]; // bottom left
        float I22_1 = img2d_xy[ x2  ][ y2  ]; // bottom right

        return (1-y_frac) * ((1-x_frac)*I11_1 + x_frac*I12_1) + (y_frac) * ((1-x_frac)*I21_1 + x_frac*I22_1);

    }

}
