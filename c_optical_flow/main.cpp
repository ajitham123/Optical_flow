#include <cv.h>
#include <highgui.h>
#include <math.h>

void test(IplImage *,IplImage *,bool);
void get_optical_flow_edges(IplImage *,IplImage *,bool);
float** grad_angle(IplImage *dx, IplImage *dy);
IplImage* edge_sobel(IplImage *, IplImage *);
long int cross_corr(IplImage *I1_gray,IplImage *I2_gray,int j, int k,int window,int x,int y,bool if_ssd);
float* remove_end(float* arr, int sizeOfArray);
int** remove_end(int** arr, int sizeOfArray);
void drawArrow(IplImage *, CvPoint p, CvPoint q, CvScalar color,int, int, int, int);

const double PI = 3.1415926;

int main()
{
    //bool im_scale = false;
    //bool if_sub_pixel = false;
    bool graphics = true;
    //int real_FOE[2] = {100,400};

    // load images:
    IplImage *I1 = cvLoadImage("indoor.jpg",CV_8S);
    IplImage *I2 = cvLoadImage("indoor_foe_180_315.jpg",CV_8S);

    get_optical_flow_edges(I1,I2,graphics);

}

void get_optical_flow_edges(IplImage *I1,IplImage *I2,bool graphics)
{
    IplImage *I1_gray = cvCreateImage(cvGetSize(I1),IPL_DEPTH_8U,1);
    IplImage *I2_gray = cvCreateImage(cvGetSize(I1),IPL_DEPTH_8U,1);
    IplImage *dx = cvCreateImage(cvGetSize(I1),IPL_DEPTH_16S,1);
    IplImage *dy = cvCreateImage(cvGetSize(I1),IPL_DEPTH_16S,1);

    cvCvtColor(I1,I1_gray,CV_RGB2GRAY);
    cvCvtColor(I2,I2_gray,CV_RGB2GRAY);
    CvSize siz = cvGetSize(I1_gray);

    cvSobel(I1_gray,dx,1,0,3);
    cvSobel(I1_gray,dy,0,1,3);
    float** Gdir = grad_angle(dx,dy);

    int kki = 250;
    for (int kkj=0;kkj<1;kkj++)
    {schar dxx = CV_IMAGE_ELEM(dx, schar,kki,2*kkj);
    schar dyy = CV_IMAGE_ELEM(dy, schar,kki,2*kkj);
    printf("dx=%d, dy=%d, angle=%f \n",dxx,dyy,Gdir[kki][kkj]);
    }

    IplImage *I1_edge = cvCreateImage(cvGetSize(I1),IPL_DEPTH_8U,1);
    I1_edge = edge_sobel(dx,dy);

    // Evaluate optical flow perpendicular to all edges
    int range = 20, tot_points =0, **point = 0;
    float *flow_mag = (float *)malloc(sizeof(float));
    float *angle    = (float *)malloc(sizeof(float));
    float flow_thresh = sqrt(2)*range;
    //printf("%d ",CV_IMAGE_ELEM(I1_edge,uchar,180,620)==255);


    for(int j=0;j<siz.height;j++)
    {
        for(int k=0;k<siz.width;k++)
        {
            //printf(" %u ",CV_IMAGE_ELEM(I1_edge,uchar,j,k));
            if(CV_IMAGE_ELEM(I1_edge,uchar,j,k)==255)   // for all edges
            {
                flow_mag = (float *)realloc(flow_mag,(tot_points+1)*sizeof(float));
                angle = (float *)realloc(angle,(tot_points+1)*sizeof(float));
                point = (int **)realloc(point,(tot_points+1)*sizeof(*point));
                point[tot_points] = (int*)malloc(2*sizeof(int));

                angle[tot_points] = Gdir[j][k]*PI/180;
                point[tot_points][0] = j;       // image coordinates
                point[tot_points][1] = k;
                //printf("\n %d, %d ",j,k);
                //printf("\n angle %f ",angle[tot_points]);
                //printf(" %d ",CV_IMAGE_ELEM(I1_edge,uchar,j,k)*CV_IMAGE_ELEM(I1_edge,uchar,j,k));


            // find search points along the gradient direction - all points
            // found will fall inside a cube of size 'range' centered about
            // 'point'
                int xt[range+1],yt[range+1],xb[range+1],yb[range+1];
                if ((angle[tot_points]>=-PI/4 && angle[tot_points]<0) || \
                        (angle[tot_points]>=3*PI/4 && angle[tot_points]<PI+0.001))
                {
                    for(int i=0;i<range+1;i++)
                    {
                        xt[i] = point[tot_points][1]-i;
                        yt[i] = round(tan(-angle[tot_points])*(float)(xt[i]-point[tot_points][1])+(float)point[tot_points][0]);
                        xb[i] = point[tot_points][1]+i;
                        yb[i] = round(tan(-angle[tot_points])*(float)(xb[i]-point[tot_points][1])+(float)point[tot_points][0]);
                    }

                }
                else if ((angle[tot_points]>=0 && angle[tot_points]<PI/4) || \
                        (angle[tot_points]>=-PI-0.001 && angle[tot_points]<-3*PI/4))
                {
                    for(int i=0;i<range+1;i++)
                    {
                        xb[i] = point[tot_points][1]-i;
                        yb[i] = round(tan(-angle[tot_points])*(float)(xb[i]-point[tot_points][1])+(float)point[tot_points][0]);
                        xt[i] = point[tot_points][1]+i;
                        yt[i] = round(tan(-angle[tot_points])*(float)(xt[i]-point[tot_points][1])+(float)point[tot_points][0]);
                    }
                }
                else
                {
                    for(int i=0;i<range+1;i++)
                    {
                        yb[i] = point[tot_points][0]+i;
                        xb[i] = round(((float)(yb[i]-point[tot_points][0])/tan(-angle[tot_points]))+(float)point[tot_points][1]);
                        yt[i] = point[tot_points][0]-i;
                        xt[i] = round(((float)(yt[i]-point[tot_points][0])/tan(-angle[tot_points]))+(float)point[tot_points][1]);
                    }
//                    printf("\n xb = ");
//                for(int i=0;i<range+1;i++)
//                {
//                    printf(" %d ",xb[i]);
//                }
                }

                // SSD Correlation along the line and evaluate optical flow

                int window = 15,best_match[2] = {0,0};
                int best_dir = 1, best_i = 0;
                long int best_corr = INFINITY, all_corr[2*range+1];
                for(int i=0;i<range+1;i++)
                {
                    // check if the search boundaries doesn't go out of the image
                    if (yb[i]-(window-1)/2>0 && yb[i]+(window-1)/2<=siz.height && \
                        xb[i]-(window-1)/2>0 && xb[i]+(window-1)/2<=siz.width && \
                        yt[i]-(window-1)/2>0 && yt[i]+(window-1)/2<=siz.height && \
                        xt[i]-(window-1)/2>0 && xt[i]+(window-1)/2<=siz.width && \
                        point[tot_points][0]-(window-1)/2>0 && point[tot_points][0]+(window-1)/2<=siz.width && \
                        point[tot_points][1]-(window-1)/2>0 && point[tot_points][1]+(window-1)/2<=siz.width)
                    {
                        long int corr_b, corr_t, corr_d;
                        bool if_ssd = true, corr_dir;
                        int corr_match[2];
                        corr_b = cross_corr(I1_gray,I2_gray,j,k,window,xb[i],yb[i],if_ssd);
                        corr_t = cross_corr(I1_gray,I2_gray,j,k,window,xt[i],yt[i],if_ssd);
                        //printf("\n corr_b %lu, corr_t %lu ",corr_b,corr_t);

//                        for(int m=0;m<window;m++)
//                            for(int n=0;n<window;n++)
//                                corr_b = corr_b + (CV_IMAGE_ELEM(I2_gray,uchar,yb[i]-(window-1)/2+m,xb[i]-(window-1)/2+n)- \
//                                    CV_IMAGE_ELEM(I1_gray,uchar,point[tot_points][0]-(window-1)/2+m,point[tot_points][1]-(window-1)/2+n))^2;
//
//                        for(int m=0;m<window;m++)
//                            for(int n=0;n<window;n++)
//                                corr_t = corr_t + (CV_IMAGE_ELEM(I2_gray,uchar,yt[i]-(window-1)/2+m,xt[i]-(window-1)/2+n)- \
//                                    CV_IMAGE_ELEM(I1_gray,uchar,point[tot_points][0]-(window-1)/2+m,point[tot_points][1]-(window-1)/2+n))^2;

                         // store all the correlation values to plot
                        all_corr[range-i] = corr_b;
                        all_corr[range+i] = corr_t;

                        // evaluate the lowest of correlation in both directions along the gradient
                        if (corr_b < corr_t)
                        {
                            corr_dir = 0;       // negative direction
                            corr_d = corr_b;
                            corr_match[0] = xb[i]; corr_match[1] = yb[i];
                        }
                        else
                        {
                            corr_dir = 1;       // flag for top portion (positive direction)
                            corr_d = corr_t;    // correlation value in the direction
                            corr_match[0] = xt[i]; corr_match[1] = yt[i];
                        }

                        // in case a better correlation is found update best matches
                        if (corr_d < best_corr)
                        {
                            best_dir   = corr_dir;
                            best_i     = i;
                            best_corr  = corr_d;
                            best_match[0] = corr_match[0];
                            best_match[1] = corr_match[1];
                        }
                    }
                }
                // sub pixel flow to refine flow magnitudes
                //printf("\n best_corr %lu ",best_corr);




                // If best correlation is in the opposite direction of gradient, update the angle
                if (best_dir == 1 && angle[tot_points]<0)
                    angle[tot_points] = angle[tot_points]+PI;
                else if (best_dir == 0 && angle[tot_points]>=0)
                    angle[tot_points] = angle[tot_points]-PI;


                flow_mag[tot_points] = (float)sqrt(pow(best_match[1]-point[tot_points][0],2)+ \
                                            pow(best_match[0]-point[tot_points][1],2));
                //printf("\n flow_mag = %f ",flow_mag[tot_points]);

                tot_points = tot_points + 1;


                // in case flow magnitude is too large - particularly near image
                // boundary - remove the flow
                if (flow_mag[tot_points-1]>flow_thresh)// || best_corr==INFINITY)
                   {
                       flow_mag = remove_end(flow_mag,tot_points);
                       angle = remove_end(angle,tot_points);
                       point = remove_end(point,tot_points);
                       tot_points = tot_points - 1;
                       //printf("\n best_corr %lu ",best_corr);
                   }



            }
        }
    }
    if(graphics)
    {
        float flow_scale = 1;
        int flow_stride = 15;
        for(int i=0;i<tot_points;i=i+flow_stride)
            {
                cvCircle(I1,cvPoint(point[i][1],point[i][0]),2,cvScalar(0,255,0),0);
                cvLine(I1,cvPoint(point[i][1],point[i][0]),cvPoint(point[i][1]+flow_scale*flow_mag[i]*cos(angle[i]), \
                                  point[i][0]-flow_scale*flow_mag[i]*sin(angle[i])),cvScalar(0,0,255),1,8,0);
                //drawArrow(I1,cvPoint(point[i][1],point[i][0]),cvPoint(point[i][1]+flow_mag[i]*cos(angle[i]), \
                                  point[i][0]-flow_mag[i]*sin(angle[i])),cvScalar(0,0,255),9,1,8,0);

                //printf("\n %d, %d",point[i][0],point[i][1]);
                printf("\n i = %d, flow = %f",i+1,flow_mag[i]);
                //cvLine(I1_gray,cvPoint(100,100),cvPoint(200,200),cvScalar(0,0,255),2,8,0);
            }
        //cvLine(I1_gray,cvPoint(100,100),cvPoint(200,200),cvScalar(0,0,255),2,8,0);
        //cvCircle(I1_gray,cvPoint(100,100),5,cvScalar(0,255,0),0);
        printf("\n\n %d",tot_points);
        cvShowImage("new", I1);
        //cvShowImage("image", angle);
        cvWaitKey(0);
    }

}
int** remove_end(int** arr, int sizeOfArray)
{
    int** temp = (int **)malloc((sizeOfArray-1)*sizeof(*arr));
    if (sizeOfArray != 0)
        memcpy(temp, arr, (sizeOfArray - 1) * sizeof(*arr));
    free(arr);
    return temp;
}

float* remove_end(float* arr, int sizeOfArray)
{
    float* temp = (float *)malloc((sizeOfArray - 1) * sizeof(float));
    if (sizeOfArray != 0)
        memcpy(temp, arr, (sizeOfArray - 1) * sizeof(float));
    free(arr);
    return temp;
}

long int cross_corr(IplImage *I1_gray,IplImage *I2_gray,int j,int k,int window,int x,int y,bool if_ssd)
{
    long int corr = 0;
    if(if_ssd)
    {
        // SSD correlation
        for(int m=0;m<window;m++)
            for(int n=0;n<window;n++)
            {
                uchar c_diff = CV_IMAGE_ELEM(I2_gray,uchar,y-(window-1)/2+m,x-(window-1)/2+n)- \
                    CV_IMAGE_ELEM(I1_gray,uchar,j-(window-1)/2+m,k-(window-1)/2+n);
                corr = corr + c_diff*c_diff;
                //printf("\nm = %d, n = %d, c_diff = %u, corr = %lu ", m,n,c_diff,corr);
            }
    }
    else
    {
        // SAD correlation
        for(int m=0;m<window;m++)
            for(int n=0;n<window;n++)
                corr = corr + abs(CV_IMAGE_ELEM(I2_gray,uchar,y-(window-1)/2+m,x-(window-1)/2+n)- \
                    CV_IMAGE_ELEM(I1_gray,uchar,j-(window-1)/2+m,k-(window-1)/2+n));

    }
   return corr;
}

float** grad_angle(IplImage *dx, IplImage *dy)
{
    // Evaluate gradient direction
    CvSize siz = cvGetSize(dx);
    float** angle = (float**)malloc(siz.height*sizeof(float*));
    for(int i=0;i<siz.height;i++)
    {
        angle[i] = (float*)malloc(siz.width*sizeof(float));
        for(int j=0;j<siz.width;j++)
        {
            schar dxx = CV_IMAGE_ELEM(dx, schar,i,2*j); // 2* because it is signed
            schar dyy = CV_IMAGE_ELEM(dy, schar,i,2*j);
            angle[i][j] = atan2(-dyy,dxx)*180/PI;
            //printf("%f  ",angle[i][j]);
        }
    }
    return angle;
}

IplImage* edge_sobel(IplImage *dx, IplImage *dy)
{
    // Edge detection using Sobel filter
    IplImage *adx = cvCreateImage(cvGetSize(dx),IPL_DEPTH_8U,1);
    IplImage *ady = cvCreateImage(cvGetSize(dx),IPL_DEPTH_8U,1);
    IplImage *ADX = cvCreateImage(cvGetSize(dx),IPL_DEPTH_8U,1);
    IplImage *ADY = cvCreateImage(cvGetSize(dx),IPL_DEPTH_8U,1);
    IplImage *I1_edge = cvCreateImage(cvGetSize(dx),IPL_DEPTH_8U,1);
    CvSize siz = cvGetSize(dx);

    cvConvertScaleAbs(dx,adx);
    cvConvertScaleAbs(dy,ady);
    CvScalar madx = cvAvg(adx), mady = cvAvg(ady);
    float edge_factor = 4;
    for(int i=0;i<siz.height;i++)
        for(int j=0;j<siz.width;j++)
        {
            ADX->imageData[i*ADX->widthStep+j*ADX->nChannels+0] = \
                CV_IMAGE_ELEM(adx, uchar,i,j) >=edge_factor*madx.val[0];
            ADY->imageData[i*ADY->widthStep+j*ADY->nChannels+0] = \
                CV_IMAGE_ELEM(ady, uchar,i,j) >=edge_factor*mady.val[0];
            I1_edge->imageData[i*I1_edge->widthStep+j*I1_edge->nChannels+0] = \
                255*(CV_IMAGE_ELEM(ADX, uchar,i,j) | CV_IMAGE_ELEM(ADY, uchar,i,j));
        }
    return I1_edge;
}

void drawArrow(IplImage *image, CvPoint p, CvPoint q, CvScalar color, int arrowMagnitude = 9, int thickness=1, int line_type=8, int shift=0)
{
    //Draw the principle line
    cvLine(image, p, q, color, thickness, line_type, shift);
    //const double PI = 3.141592653;
    //compute the angle alpha
    double angle = atan2((double)p.y-q.y, (double)p.x-q.x);
    //compute the coordinates of the first segment
    p.x = (int) ( q.x +  arrowMagnitude * cos(angle + PI/4));
    p.y = (int) ( q.y +  arrowMagnitude * sin(angle + PI/4));
    //Draw the first segment
    cvLine(image, p, q, color, thickness, line_type, shift);
    //compute the coordinates of the second segment
    p.x = (int) ( q.x +  arrowMagnitude * cos(angle - PI/4));
    p.y = (int) ( q.y +  arrowMagnitude * sin(angle - PI/4));
    //Draw the second segment
    cvLine(image, p, q, color, thickness, line_type, shift);
}

void test(IplImage *I1,IplImage *I2,bool graphics)
{

    //cvSmooth( I1, I1, CV_GAUSSIAN, 3,0,0,0);
    IplImage *I1_gray = cvCreateImage(cvGetSize(I1),IPL_DEPTH_8U,1);
    cvCvtColor(I1,I1_gray,CV_RGB2GRAY);
    CvSize siz = cvGetSize(I1_gray);
    //printf("Width %d and Height %d",siz.width,siz.height);
    IplImage *dx = cvCreateImage(cvGetSize(I1_gray),IPL_DEPTH_8U,1);
    IplImage *dy = cvCreateImage(cvGetSize(I1_gray),IPL_DEPTH_8U,1);

    float kernel_x[9] = {0,0,0,-1.0,0,1.0,0,0,0};
    float kernel_y[9] = {0,-1.0,0,0,0,0,0,1.0,0};
    //float kernel_x[9] = {0.25/4,0.5/4,0.25/4,.5/4,1/4,.5/4,0.25/4,0.5/4,0.25/4};    // gaussian kernel
    CvMat filter_x,filter_y;
    //kernel_y cvMat(3,3,CV_32F, {0,-1,0,0,0,0,0,1,0});

    filter_x = cvMat(3,3,CV_32F, kernel_x);
    filter_y = cvMat(3,3,CV_32F, kernel_y);

    cvFilter2D(I1_gray,dx,&filter_x);
    cvFilter2D(I1_gray,dy,&filter_y);

    cvNamedWindow("image");
    cvShowImage("image", dx);
    cvWaitKey(0);
    cvReleaseImage(&dx);
    cvDestroyAllWindows();
}
