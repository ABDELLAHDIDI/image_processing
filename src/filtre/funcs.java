/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package classes;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import org.opencv.core.Mat;
import org.opencv.core.Core;
import org.opencv.core.Scalar;
import org.opencv.core.Size;
import org.opencv.imgproc.Imgproc;

/**
 *
 * @author SohaibID
 */
public class funcs {
    
    public static double[][] getMkernel(int size){
        double[][] res= new double[size][size];
        int sum=size*size;
        for(int i=0;i<size;i++){
            for(int j=0;j<size;j++){
                res[i][j]=1.0/sum;
            }
        }
        return res;
    }
    
    
    public static double[][] getGaussianKernel(int size, double sigma) {
        double[][] kernel = new double[size][size];
        double sum = 0;

        int radius = size / 2;
        for (int y = -radius; y <= radius; y++) {
            for (int x = -radius; x <= radius; x++) {
                double value = Math.exp(-(x * x + y * y) / (2 * sigma * sigma));
                kernel[y + radius][x + radius] = value;
                sum += value;
            }
        }

        // Normalize the kernel
        for (int y = 0; y < size; y++) {
            for (int x = 0; x < size; x++) {
                kernel[y][x] /= sum;
            }
        }

        return kernel;
    }
    
    public static Mat getPaddedImage(Mat inputImage ,int paddingSize){
        var rows_nbr=inputImage.rows();
        var cols_nbr=inputImage.cols();
        
        var paddedImage = new Mat(rows_nbr + 2*paddingSize, cols_nbr + 2*paddingSize, inputImage.type());
        paddedImage.setTo(new Scalar(0));
        
        var sub=paddedImage.submat(paddingSize, paddedImage.rows() - paddingSize, paddingSize, paddedImage.cols() - paddingSize);
        inputImage.copyTo(sub);
//        for (int i = 0; i < rows_nbr; i++) {
//            for (int j = 0; j < cols_nbr; j++) {
//                paddedImage.put(i + paddingSize, j + paddingSize, inputImage.get(i, j));
//            }
//        }
        return paddedImage;
    }
    
    public static Mat toGray(Mat inputImage){
        Mat outputImage = new Mat(inputImage.size(), inputImage.type(),new Scalar(0));
        var rows_nbr=outputImage.rows();
        var cols_nbr=outputImage.cols();
        for (int i = 0; i < rows_nbr; i++) {
            for (int j = 0; j < cols_nbr; j++) {
                var imagePixel=inputImage.get(i, j);
                var mean=imagePixel[0]+imagePixel[1]+imagePixel[2];
                mean=mean/3;
                imagePixel[0]=mean;
                imagePixel[1]=mean;
                imagePixel[2]=mean;
                outputImage.put(i, j, imagePixel);
            }
        }
        
        return outputImage;
        
    }
    
    public static Mat binary(Mat inputImage){
        Mat outputImage = new Mat(inputImage.size(), inputImage.type(),new Scalar(0));
        var rows_nbr=outputImage.rows();
        var cols_nbr=outputImage.cols();
        var chanels=outputImage.channels();
        for (int i = 0; i < rows_nbr; i++) {
            for (int j = 0; j < cols_nbr; j++) {
                var imagePixel=inputImage.get(i, j);
                var mean=imagePixel[0]+imagePixel[1]+imagePixel[2];
                mean=mean/3;
                if(mean>150)mean=255;
                else mean=0;
                imagePixel[0]=mean;
                imagePixel[1]=mean;
                imagePixel[2]=mean;
                if(chanels==4){
                    imagePixel[3]=255;
                }
                outputImage.put(i, j, imagePixel);
            }
        }
        
        return outputImage;
        
    }
    
    public static void filterImage(Mat inputImage ,Mat outputImage,double [][] kernelMat){
        var kernelSize=kernelMat.length;
        var marge=kernelSize/2;
        var rows_nbr=outputImage.rows();
        var cols_nbr=outputImage.cols();
        var chanels=outputImage.channels();
        double[] pixel0=outputImage.get(0,0).clone();
        
        var inputmat=new double[rows_nbr][cols_nbr][chanels];
        for (int i = 0; i < rows_nbr; i++) {
            for (int j = 0; j < cols_nbr; j++) {
                inputmat[i][j]=inputImage.get(i, j).clone();
                
            }
        }
        
        
        for (int i = 0; i < rows_nbr; i++) {
            for (int j = 0; j < cols_nbr; j++) {
                var sum=pixel0.clone();
                int count = 0;
                for (int k = 0; k < kernelSize; k++) {
                    int row = i + k - marge;
                    if(row<0 || row>=rows_nbr)continue;
                    for (int l = 0; l < kernelSize; l++) {
                        
                        int col = j + l - marge;
                        if(col>=0 && col<cols_nbr) {
                            var imagePixel=inputmat[row][col];
                            var kernelValue=kernelMat[k][l];
                            count++;
                            for(int e=0; e<chanels ; e++){
                                sum[e]+=kernelValue*imagePixel[e];
                                
                            }
                        }
                    }
                }
                
                
                outputImage.put(i, j, sum);
            }
        }
    }
    
    
    
    public static void filterImageWP(Mat inputImage ,Mat outputImage,double[][] kernelMat){
        var kernelSize=kernelMat.length;
        int paddingSize = kernelSize/2;
        int count=kernelSize*kernelSize;
        // get the padded image
        Mat paddedImage = getPaddedImage(inputImage, paddingSize);
        var rows_nbr=outputImage.rows();
        var cols_nbr=outputImage.cols();
        var chanels=outputImage.channels();
        
        var prows_nbr=paddedImage.rows();
        var pcols_nbr=paddedImage.cols();
        var inputmat=new double[prows_nbr][pcols_nbr][chanels];
        for (int i = 0; i < prows_nbr; i++) {
            for (int j = 0; j < pcols_nbr; j++) {
                inputmat[i][j]=paddedImage.get(i, j).clone();
            }
        }
        
        double[] pixel0=outputImage.get(0,0).clone();
        for (int i = 0; i < rows_nbr; i++) {
            for (int j = 0; j < cols_nbr; j++) {
                var sum=pixel0.clone();
                for (int k = 0; k < kernelSize; k++) {
                    for (int l = 0; l < kernelSize; l++) {
                        var kernelValue=kernelMat[k][l];     
                        var imagePixel=inputmat[i + k][j + l];
                        for(int e=0; e<chanels ; e++){
                            sum[e]+=kernelValue*imagePixel[e];
                        }
                        
                    }
                }
                
                
                outputImage.put(i, j, sum);
            }
        }
    }
    
    public static BufferedImage medianBlur(BufferedImage inputImage, int kernelSize) {
          int width = inputImage.getWidth();
        int height = inputImage.getHeight();

        BufferedImage outputImage = new BufferedImage(width, height, inputImage.getType());

        int paddingSize = (kernelSize - 1) / 2;
        
        int[][] inputMat=new int[width][height];
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                inputMat[x][y]=inputImage.getRGB(x, y);
            }
        }
        
        // Apply median blur
        for (int y = paddingSize; y < height - paddingSize; y++) {
            for (int x = paddingSize; x < width - paddingSize; x++) {
                int[] neighbors = new int[kernelSize * kernelSize];
                int index = 0;

                // Extract neighborhood pixel values
                for (int ky = y - paddingSize; ky <= y + paddingSize; ky++) {
                    for (int kx = x - paddingSize; kx <= x + paddingSize; kx++) {
                        Color pixelColor = new Color(inputMat[kx][ky]);
                        neighbors[index++] = pixelColor.getRGB();
                    }
                }

                // Sort the neighborhood pixel values
                Arrays.sort(neighbors);

                // Calculate the median
                int median;
                if (kernelSize % 2 == 0) {
                    median = (neighbors[kernelSize * kernelSize / 2 - 1] + neighbors[kernelSize * kernelSize / 2]) / 2;
                } else {
                    median = neighbors[kernelSize * kernelSize / 2];
                }

                outputImage.setRGB(x, y, median);
            }
        }

        return outputImage;
    }
    
    public static void filterM(Mat inputImage ,Mat outputImage,int kernelSize){
//        Imgproc.medianBlur(inputImage, outputImage, kernelSize);
        Imgproc.blur(inputImage, outputImage, new Size(kernelSize,kernelSize));
    }
    
    
  public static void pixelate(Mat inputImage ,Mat outputImage,int kernelSize){
        
        var marge=kernelSize/2;
        var rows_nbr=outputImage.rows();
        var cols_nbr=outputImage.cols();
        
        for (int i = marge; i < rows_nbr; i+=2*marge) {
            for (int j = marge; j < cols_nbr; j+=2*marge) {
                var pixel=inputImage.get(i, j);
                var rowStart=i-marge;
                var rowEnd=i+marge;
                var colStart=j-marge;
                var colEnd=j+marge;
                if(rowEnd>rows_nbr)rowEnd=rows_nbr;
                if(colEnd>cols_nbr)colEnd=cols_nbr;
                var subM=outputImage.submat(rowStart, rowEnd, colStart, colEnd);
                subM.setTo(new Scalar(pixel));
            }
        }
        
    }
    
    
}
