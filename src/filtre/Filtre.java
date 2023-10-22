/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Main.java to edit this template
 */
package filtre;

import java.awt.Color;
import java.awt.Image;
import java.awt.Point;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Stack;
import javax.imageio.ImageIO;

/**
 *
 * @author ABDELLAH DIDI
 */
public class Filtre {
    
    
    
        public BufferedImage   Filtre_sobel (Image img,int  kernelSize ) {
 
   
        double [][] sobelKernelX = generateSobelKernelX(kernelSize);
        double [][] sobelKernelY = generateSobelKernelY(kernelSize);

        
            // Charger l'image à filtrer
            BufferedImage image = (BufferedImage) img;

            // Convertir l'image en niveaux de gris
            BufferedImage grayImage = convertToGray(image);

            // Appliquer le filtre de Sobel pour détecter les contours horizontaux et verticaux
   BufferedImage sobelImage = applySobelFilter(grayImage, sobelKernelX, sobelKernelY);
   return sobelImage;
    }
        
    private static double [][] generateSobelKernelX(int size) {
        double [][] kernelX = new double [size][size];
        double  middle = size / 2;

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (j == middle) {
                    kernelX[i][j] = 2 * (i - middle);
                } else {
                    kernelX[i][j] = (i - middle) * (j - middle);
                }
            }
        }

        return kernelX;
    }

    private static double [][] generateSobelKernelY(int size) {
        double [][] kernelY = new double [size][size];
        double  middle = size / 2;

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (i == middle) {
                    kernelY[i][j] = 2 * (j - middle);
                } else {
                    kernelY[i][j] = (i - middle) * (j - middle);
                }
            }
        }

        return kernelY;
    }
    
    private static BufferedImage convertToGray(BufferedImage image) {
        BufferedImage grayImage = new BufferedImage(image.getWidth(), image.getHeight(), BufferedImage.TYPE_BYTE_GRAY);
        grayImage.getGraphics().drawImage(image, 0, 0, null);
        return grayImage;
    }

    private static BufferedImage applySobelFilter(BufferedImage image, double [][] sobelXKernel, double [][] sobelYKernel) {
        int width = image.getWidth();
        int height = image.getHeight();
        BufferedImage filteredImage = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY);

        int kernelSize = sobelXKernel.length;
        int kernelOffset = kernelSize / 2;

        for (int x = kernelOffset; x < width - kernelOffset; x++) {
            for (int y = kernelOffset; y < height - kernelOffset; y++) {
                int gx = 0;
                int gy = 0;

                for (int i = 0; i < kernelSize; i++) {
                    for (int j = 0; j < kernelSize; j++) {
                        int pixelValue = new Color(image.getRGB(x + i - kernelOffset, y + j - kernelOffset)).getRed();

                        gx += sobelXKernel[i][j] * pixelValue;
                        gy += sobelYKernel[i][j] * pixelValue;
                    }
                }

                int magnitude = (int) Math.sqrt(gx * gx + gy * gy);
                magnitude = Math.min(magnitude, 255);  // Limiter les valeurs à 255 pour éviter les débordements

                int grayValue = new Color(magnitude, magnitude, magnitude).getRGB();
                filteredImage.setRGB(x, y, grayValue);
            }
        }

        return filteredImage;
    }
    

    
    
    
         public BufferedImage   Filtre_Canny (Image img,int  kernelSize ,int Gaussiansize) {
 
   
        double [][] sobelKernelX = generateSobelKernelX(kernelSize);
        double [][] sobelKernelY = generateSobelKernelY(kernelSize);

        
               BufferedImage image = (BufferedImage) img;

            // Convertir l'image en niveaux de gris
            BufferedImage grayImage = convertToGray(image);

            // Appliquer le filtre de Canny
   BufferedImage cannyImage = applyCannyFilter(grayImage, sobelKernelX, sobelKernelY,Gaussiansize);
   
   return cannyImage;
    }
    
    private static BufferedImage applyCannyFilter(BufferedImage image, double [][] sobelKernelX, double [][] sobelKernelY,int size) {
        int width = image.getWidth();
        int height = image.getHeight();
        BufferedImage cannyImage = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY);

        // Définir le noyau pour le lissage (gaussian kernel)
   
        
        
       


double  [][] gaussianKernel = generateGaussianKernel(size, 1);


        // Appliquer le lissage
        BufferedImage smoothedImage = applyConvolution(image, gaussianKernel);

     
        // Calculer les gradients horizontaux et verticaux
        BufferedImage gradientXImage = applyConvolution(smoothedImage, sobelKernelX);
        BufferedImage gradientYImage = applyConvolution(smoothedImage, sobelKernelY);

        // Calculer la magnitude du gradient
        BufferedImage magnitudeImage = calculateGradientMagnitude(gradientXImage, gradientYImage);

        // Appliquer la suppression des non-maximaux
        BufferedImage suppressedImage = suppressNonMaxima(magnitudeImage, gradientXImage, gradientYImage);

        // Appliquer la détection des contours par seuillage
        BufferedImage thresholdedImage = applyThresholding(suppressedImage, 30, 100);

        // Appliquer la hysteresis pour relier les contours
        BufferedImage cannyEdgesImage = applyHysteresis(thresholdedImage, 20, 60);

        return cannyEdgesImage;
    }

    private static BufferedImage applyConvolution(BufferedImage image, double [][] kernel) {
        int width = image.getWidth();
        int height = image.getHeight();
        int kernelSize = kernel.length;
        int kernelHalfSize = kernelSize / 2;
        BufferedImage result = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY);

        for (int y = kernelHalfSize; y < height - kernelHalfSize; y++) {
            for (int x = kernelHalfSize; x < width - kernelHalfSize; x++) {
                int sum = 0;

                for (int i = 0; i < kernelSize; i++) {
                    for (int j = 0; j < kernelSize; j++) {
                        int pixel = new Color(image.getRGB(x + j - kernelHalfSize, y + i - kernelHalfSize)).getRed();
                        sum += kernel[i][j] * pixel;
                    }
                }

                int newValue = Math.min(Math.max(sum, 0), 255);
                Color newColor = new Color(newValue, newValue, newValue);
                result.setRGB(x, y, newColor.getRGB());
            }
        }

        return result;
    }

    private static BufferedImage calculateGradientMagnitude(BufferedImage gradientXImage, BufferedImage gradientYImage) {
        int width = gradientXImage.getWidth();
        int height = gradientXImage.getHeight();
        BufferedImage magnitudeImage = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY);

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int gradientX = new Color(gradientXImage.getRGB(x, y)).getRed();
                int gradientY = new Color(gradientYImage.getRGB(x, y)).getRed();
                int magnitude = (int) Math.sqrt(gradientX * gradientX + gradientY * gradientY);
                int clampedMagnitude = Math.min(magnitude, 255);
                Color magnitudeColor = new Color(clampedMagnitude, clampedMagnitude, clampedMagnitude);

//                Color magnitudeColor = new Color(magnitude, magnitude, magnitude);
                magnitudeImage.setRGB(x, y, magnitudeColor.getRGB());
            }
        }

        return magnitudeImage;
    }

    private static BufferedImage suppressNonMaxima(BufferedImage magnitudeImage, BufferedImage gradientXImage, BufferedImage gradientYImage) {
        int width = magnitudeImage.getWidth();
        int height = magnitudeImage.getHeight();
        BufferedImage suppressedImage = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY);

        for (int y = 1; y < height - 1; y++) {
            for (int x = 1; x < width - 1; x++) {
                int magnitude = new Color(magnitudeImage.getRGB(x, y)).getRed();
                int gradientX = new Color(gradientXImage.getRGB(x, y)).getRed();
                int gradientY = new Color(gradientYImage.getRGB(x, y)).getRed();

                int orientation = (int) (Math.atan2(gradientY, gradientX) * (180 / Math.PI));
                if (orientation < 0) {
                    orientation += 180;
                }

                int neighbor1 = 0;
                int neighbor2 = 0;

                if ((orientation >= 0 && orientation < 22.5) || (orientation >= 157.5 && orientation <= 180)) {
                    neighbor1 = new Color(magnitudeImage.getRGB(x + 1, y)).getRed();
                    neighbor2 = new Color(magnitudeImage.getRGB(x - 1, y)).getRed();
                } else if (orientation >= 22.5 && orientation < 67.5) {
                    neighbor1 = new Color(magnitudeImage.getRGB(x + 1, y - 1)).getRed();
                    neighbor2 = new Color(magnitudeImage.getRGB(x - 1, y + 1)).getRed();
                } else if (orientation >= 67.5 && orientation < 112.5) {
                    neighbor1 = new Color(magnitudeImage.getRGB(x, y - 1)).getRed();
                    neighbor2 = new Color(magnitudeImage.getRGB(x, y + 1)).getRed();
                } else if (orientation >= 112.5 && orientation < 157.5) {
                    neighbor1 = new Color(magnitudeImage.getRGB(x - 1, y - 1)).getRed();
                    neighbor2 = new Color(magnitudeImage.getRGB(x + 1, y + 1)).getRed();
                }

                if (magnitude >= neighbor1 && magnitude >= neighbor2) {
                    suppressedImage.setRGB(x, y, new Color(magnitude, magnitude, magnitude).getRGB());
                } else {
                    suppressedImage.setRGB(x, y, new Color(0, 0, 0).getRGB());
                }
            }
        }

        return suppressedImage;
    }

    private static BufferedImage applyThresholding(BufferedImage image, int lowThreshold, int highThreshold) {
        int width = image.getWidth();
        int height = image.getHeight();
        BufferedImage thresholdedImage = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY);

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int magnitude = new Color(image.getRGB(x, y)).getRed();

                if (magnitude >= highThreshold) {
                    thresholdedImage.setRGB(x, y, new Color(255, 255, 255).getRGB());
                } else if (magnitude >= lowThreshold) {
                    thresholdedImage.setRGB(x, y, new Color(128, 128, 128).getRGB());
                } else {
                    thresholdedImage.setRGB(x, y, new Color(0, 0, 0).getRGB());
                }
            }
        }

        return thresholdedImage;
    }

    private static BufferedImage applyHysteresis(BufferedImage image, int lowThreshold, int highThreshold) {
        int width = image.getWidth();
        int height = image.getHeight();
        BufferedImage hysteresisImage = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY);
        boolean[][] visited = new boolean[width][height];

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int magnitude = new Color(image.getRGB(x, y)).getRed();

                if (magnitude >= highThreshold && !visited[x][y]) {
                   
                   traceEdges(image, hysteresisImage, visited, x, y, lowThreshold);
           
                // traceEdges(image, hysteresisImage, visited,  lowThreshold);
                }
            }
        }

        return hysteresisImage;
    }

    private static void traceEdges(BufferedImage image, BufferedImage hysteresisImage, boolean[][] visited, int x, int y, int threshold) {
        int width = image.getWidth();
        int height = image.getHeight();

        if (x >= 0 && x < width && y >= 0 && y < height && !visited[x][y]) {
            visited[x][y] = true;
            int magnitude = new Color(image.getRGB(x, y)).getRed();

            if (magnitude >= threshold) {
                hysteresisImage.setRGB(x, y, new Color(255, 255, 255).getRGB());

                
                    traceEdges(image, hysteresisImage, visited, x - 1, y - 1, threshold);
                traceEdges(image, hysteresisImage, visited, x, y - 1, threshold);
                traceEdges(image, hysteresisImage, visited, x + 1, y - 1, threshold);
                traceEdges(image, hysteresisImage, visited, x - 1, y, threshold);
                traceEdges(image, hysteresisImage, visited, x + 1, y, threshold);
                traceEdges(image, hysteresisImage, visited, x - 1, y + 1, threshold);
                traceEdges(image, hysteresisImage, visited, x, y + 1, threshold);
                traceEdges(image, hysteresisImage, visited, x + 1, y + 1, threshold);
           
            }
        }
    }

    public static double [][] generateGaussianKernel(int size, int  sigma) {
    double  [][] kernel = new double  [size][size];
    double   sum = 0.0;

    int halfSize = size / 2;

    for (int y = -halfSize; y <= halfSize; y++) {
        for (int x = -halfSize; x <= halfSize; x++) {
            double exponent = -(x * x + y * y) / (2 * sigma * sigma);
            double value = Math.exp(exponent) / (2 * Math.PI * sigma * sigma);
            kernel[y + halfSize][x + halfSize] = value;
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







        public BufferedImage   Filtre_Prewitt (Image img,int  Masksize ) {
            // Charger l'image à filtrer
            BufferedImage image = (BufferedImage) img;


            // Appliquer le filtre de Sobel pour détecter les contours horizontaux et verticaux
   BufferedImage PrewittImage = applyFilter(image,Masksize);
   return PrewittImage;
    }
  
  public static BufferedImage applyFilter(BufferedImage image,int size) {
        int width = image.getWidth();
        int height = image.getHeight();

        BufferedImage resultImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);

        int[][] horizontalMask = generatePrewittHorizontalMask( size);
        int[][] verticalMask =  generatePrewittVerticalMask( size);
        
        for (int y = 1; y < height - 1; y++) {
            for (int x = 1; x < width - 1; x++) {
                int horizontalGradient = applyMask(image, x, y, horizontalMask);
                int verticalGradient = applyMask(image, x, y, verticalMask);
                int gradientMagnitude = (int) Math.sqrt(horizontalGradient * horizontalGradient + verticalGradient * verticalGradient);

//                resultImage.setRGB(x, y, new Color(gradientMagnitude, gradientMagnitude, gradientMagnitude).getRGB());
int normalizedMagnitude = Math.min(255, Math.max(0, gradientMagnitude));
resultImage.setRGB(x, y, new Color(normalizedMagnitude, normalizedMagnitude, normalizedMagnitude).getRGB());

            }
        }

        return resultImage;
    }

  private static int applyMask(BufferedImage image, int x, int y, int[][] mask) {
        int result = 0;

        for (int i = -1; i <= 1; i++) {
            for (int j = -1; j <= 1; j++) {
                int pixelValue = new Color(image.getRGB(x + j, y + i)).getRed();
                result += pixelValue * mask[i + 1][j + 1];
            }
        }

        return result;
    }

  private static int[][] generatePrewittHorizontalMask(int size) {
    int[][] mask = new int[size][size];
    int center = size / 2;

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (j == center) {
                mask[i][j] = 0;
            } else if (j < center) {
                mask[i][j] = -1;
            } else {
                mask[i][j] = 1;
            }
        }
    }

    return mask;
}

  private static int[][] generatePrewittVerticalMask(int size) {
    int[][] mask = new int[size][size];
    int center = size / 2;

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (i == center) {
                mask[i][j] = 0;
            } else if (i < center) {
                mask[i][j] = -1;
            } else {
                mask[i][j] = 1;
            }
        }
    }

    return mask;
}
  
  
  
  
     public static BufferedImage   Filtre_Laplacian (Image img,int  size ) {
      
         
            BufferedImage outputImage = applyLaplacianFilter((BufferedImage) img ,size);
         
   return outputImage;
    }
  
    private static BufferedImage applyLaplacianFilter(BufferedImage img,int size) {
         int[][] laplacianKernel = createLaplacianKernel(size);
           
        
        // Load the input image
        BufferedImage inputImage =(BufferedImage)img ;

        // Create a grayscale version of the image
        BufferedImage grayImage = new BufferedImage(
            inputImage.getWidth(), inputImage.getHeight(), BufferedImage.TYPE_BYTE_GRAY);
        for (int y = 0; y < inputImage.getHeight(); y++) {
            for (int x = 0; x < inputImage.getWidth(); x++) {
                Color color = new Color(inputImage.getRGB(x, y));
                int grayValue = (int) (0.2126 * color.getRed() +
                                       0.7152 * color.getGreen() +
                                       0.0722 * color.getBlue());
                grayImage.setRGB(x, y, new Color(grayValue, grayValue, grayValue).getRGB());
            }
        }

        // Apply the Laplacian filter
        BufferedImage outputImage = new BufferedImage(
            grayImage.getWidth(), grayImage.getHeight(), BufferedImage.TYPE_BYTE_GRAY);
//        int[][] laplacianKernel = {
//            {-1, -1, -1},
//            {-1,  8, -1},
//            {-1, -1, -1}
//        };
        for (int y = 1; y < grayImage.getHeight() - 1; y++) {
            for (int x = 1; x < grayImage.getWidth() - 1; x++) {
                int sum = 0;
                for (int ky = -1; ky <= 1; ky++) {
                    for (int kx = -1; kx <= 1; kx++) {
                        int pixelValue = new Color(grayImage.getRGB(x + kx, y + ky)).getRed();
                        sum += pixelValue * laplacianKernel[ky + 1][kx + 1];
                    }
                }
                sum = Math.min(Math.max(sum, 0), 255);
                outputImage.setRGB(x, y, new Color(sum, sum, sum).getRGB());
            }
        }

       return outputImage;
        
        
    }
    
    private static  int[][] createLaplacianKernel(int size) {
    int[][] kernel = new int[size][size];

    int center = size / 2;
    
   
        kernel[center][center] = -4;
        kernel[center-1][center] = kernel[center+1][center] = kernel[center][center-1] = kernel[center][center+1] = 1;


    return kernel;
}



}