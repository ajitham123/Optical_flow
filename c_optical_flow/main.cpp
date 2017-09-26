/*
    Author: Ajith Anil Meera, September 2017, MAV LAb, TU Delft
    C program to find the optical flow on the edge points in the direction perpendicular
    to the edges. Evaluate the Focus of Expansion from the flow by removing the outliers using RANSAC.
*/

#include <cv.h>
#include <highgui.h>
#include <math.h>
#include <stdlib.h>

/*
    Structure that stores the edge points with flow magnitudes along the angle
*/
struct opt_flow
{
    int **point, tot_points;
    float *flow_mag,*angle;
};

void test(IplImage *,IplImage *,bool);
opt_flow get_optical_flow_edges(IplImage *,IplImage *,bool);
float* find_FOE(opt_flow flow);
opt_flow find_no_flow(opt_flow flow);
float** grad_angle(IplImage *dx, IplImage *dy);
IplImage* edge_sobel(IplImage *, IplImage *);
unsigned long int cross_corr(IplImage *I1_gray,IplImage *I2_gray,int j, int k,int window,int x,int y,bool if_ssd);
float* remove_end(float* arr, int sizeOfArray);
int** remove_end(int** arr, int sizeOfArray);
float** mat_multiply(float**,int,int,float**,int,int);
float** pinv(float**,int,int);
float Determinant(float **a,int n);
float** CoFactor(float **a,int n);
float** Transpose(float **a,int m,int n);
float** mat_inv(float **a,int m);
void drawArrow(IplImage *, CvPoint p, CvPoint q, CvScalar color,int, int, int, int);
float* least_square(opt_flow flow);
float* dist2lines(opt_flow flow, float, float);
float mean(float* a, int n);
float std_dev(float* a, int n);
opt_flow RANSAC(opt_flow flow,float* error_distances);
opt_flow remove_outlier_flow(opt_flow flow, float* error_distances,float mean_err,float err_std_thresh);
void display_flow(opt_flow flow, float *);
opt_flow random_sample(opt_flow flow, int n);

const double PI = 3.1415926;


int main()
{
    //bool im_scale = false;
    //bool if_sub_pixel = false;
    bool graphics = false;
    opt_flow flow;
    //int real_FOE[2] = {100,400};

    // load images:
    IplImage *I1 = cvLoadImage("indoor.jpg",CV_8S);
    IplImage *I2 = cvLoadImage("indoor_foe_200_300.jpg",CV_8S);

    flow = get_optical_flow_edges(I1,I2,graphics);
    printf("\n Number of points %d",flow.tot_points);

    float *FOE = find_FOE(flow);
    printf("\n FOE in main = %f,%f",FOE[0],FOE[1]);

    return 0;
}

/*
    Function to find the least square solution for the intersection of lines that are
    tangent to the edge points with no flow in the direction normal to the edge.
*/
float* least_square(opt_flow flow)
{
    float **A = (float **)malloc(flow.tot_points*sizeof(float *));
    float **B = (float **)malloc(flow.tot_points*sizeof(float *));
    for (int i=0; i<flow.tot_points; i++)
    {
        A[i] = (float *)malloc(2 * sizeof(float));
        B[i] = (float *)malloc(1 * sizeof(float));

        A[i][0] = -(tan(flow.angle[i]));   A[i][1] = 1;
        B[i][0] = (float)flow.point[i][1] - (tan(flow.angle[i]))*flow.point[i][0];
    }

    float **FOE = mat_multiply(pinv(A,flow.tot_points,2),2,flow.tot_points,B,flow.tot_points,1);

    // Convert 2D array to a 1D array and return
    float *FOE1 = (float *)malloc(2*sizeof(float));
    FOE1[0] = FOE[0][0];    FOE1[1] = FOE[1][0];
    //printf("\n FOE in LS = %f, %f",FOE1[0],FOE1[1]);
    return FOE1;
}

/*
    Function to find the FOE from the flow
*/
float* find_FOE(opt_flow flow)
{
    flow = find_no_flow(flow);
    flow.flow_mag = (float *)realloc(flow.flow_mag,(flow.tot_points)*sizeof(float));
    flow.angle = (float *)realloc(flow.angle,(flow.tot_points)*sizeof(float));
    flow.point = (int **)realloc(flow.point,(flow.tot_points)*sizeof(*flow.point));
    flow.point[flow.tot_points] = (int*)malloc(2*sizeof(int));

//    flow.flow_mag = (float *)malloc(sizeof(float));
//    flow.angle =    (float *)malloc(sizeof(float));
//    flow.point = 0;
//
//    for(int i=0;i<kk_flow.tot_points;i++)
//    {
//        flow.flow_mag = (float *)realloc(flow.flow_mag,(i+1)*sizeof(float));
//        flow.angle = (float *)realloc(flow.angle,(i+1)*sizeof(float));
//        flow.point = (int **)realloc(flow.point,(i+1)*sizeof(*flow.point));
//        flow.point[i] = (int*)malloc(2*sizeof(int));
//
//        flow.flow_mag[i] = kk_flow.flow_mag[i];
//        flow.angle[i] = kk_flow.angle[i];
//        flow.point[i][0] = kk_flow.point[i][0];
//        flow.point[i][1] = kk_flow.point[i][1];
//    }
//    flow.tot_points = kk_flow.tot_points;
    //printf("\n Angle %f ",flow.angle[flow.tot_points-1]);

    //rr(flow);

    printf("\n Number of no flow vectors: %d ",flow.tot_points);



    float *FOE = least_square(flow);
    float *error_distances = dist2lines(flow,FOE[0],FOE[1]);
    //printf("\n FOE = %f, %f",FOE[0], FOE[1]);
    //printf("\n dist : %f ",error_distances[flow.tot_points-1]);
    //printf("\n mean  = %f",mean(error_distances,flow.tot_points));
    //printf("\n std   = %f",std_dev(error_distances,flow.tot_points));
    flow = RANSAC(flow,error_distances);
    FOE = least_square(flow);
    display_flow(flow,FOE);

    return FOE;
}

/*
    Function to display the results of the algorithm (flow, FOE, RANSAC steps)
*/
void display_flow(opt_flow flow, float *FOE)
{
    IplImage *I1 = cvLoadImage("indoor.jpg",CV_8S);
    for(int i=0;i<flow.tot_points;i=i+1)
        cvLine(I1,cvPoint(flow.point[i][1]+400*sin(flow.angle[i]),flow.point[i][0]+400*cos(flow.angle[i])), \
                  cvPoint(flow.point[i][1]-400*sin(flow.angle[i]),flow.point[i][0]-400*cos(flow.angle[i])), \
                          cvScalar(0,0,255),1,8,0);
    for(int i=0;i<flow.tot_points;i=i+1)
            {
                //printf("\n angles %f ",flow.angle[i]);

                cvCircle(I1,cvPoint(flow.point[i][1],flow.point[i][0]),2,cvScalar(0,255,0),0);
                //cvLine(I1,cvPoint(flow.point[i][1],flow.point[i][0]),cvPoint(flow.point[i][1]+flow.flow_mag[i]*cos(flow.angle[i]), \
                                  flow.point[i][0]-flow.flow_mag[i]*sin(flow.angle[i])),cvScalar(0,255,0),1,8,0);
                //drawArrow(I1,cvPoint(point[i][1],point[i][0]),cvPoint(point[i][1]+flow_mag[i]*cos(angle[i]), \
                                  point[i][0]-flow_mag[i]*sin(angle[i])),cvScalar(0,0,255),9,1,8,0);

                //printf("\n %d, %d",point[i][0],point[i][1]);
                //printf("\n i = %d, flow = %f",i+1,flow_mag[i]);
                //cvLine(I1_gray,cvPoint(100,100),cvPoint(200,200),cvScalar(0,0,255),2,8,0);
            }
    cvCircle(I1,cvPoint(FOE[1],FOE[0]),2,cvScalar(255,0,0),0);

    cvShowImage("image",I1);
    cvWaitKey(0);
}

/*
    Function to perform RANSAC over the set of lines so as to remove the outlier lines.
    Sum of perpendicular distances from the FOE to all the inlier lines are minimized.
    All lines farther than a threshold standard deviation are removed iteratively.
*/
opt_flow RANSAC(opt_flow flow,float* error_distances)
{
    int kk = 0;
    float err_std_thresh = 1.8;
    float mean_err = mean(error_distances,flow.tot_points);
    while(mean_err>3 && kk<70)
    {
        mean_err = mean(error_distances,flow.tot_points);
        //printf("\n Mean err : %f ",mean_err);
        //printf("\n Tot_points : %d ",flow.tot_points);
        opt_flow kk_flow  = remove_outlier_flow(flow,error_distances,mean_err,err_std_thresh);

        /*
            Updating the flow at each iteration using a temporary kk_flow object.
            Necessary to get around the dynamic memory allocation issues
        */
        flow.flow_mag = (float *)malloc(sizeof(float));
        flow.angle =    (float *)malloc(sizeof(float));
        flow.point = 0;

        // Initial implementation of RANSAC
        for(int i=0;i<kk_flow.tot_points;i++)
        {
                flow.flow_mag = (float *)realloc(flow.flow_mag,(i+1)*sizeof(float));
                flow.angle = (float *)realloc(flow.angle,(i+1)*sizeof(float));
                flow.point = (int **)realloc(flow.point,(i+1)*sizeof(*flow.point));
                flow.point[i] = (int*)malloc(2*sizeof(int));

                flow.flow_mag[i] = kk_flow.flow_mag[i];
                flow.angle[i] = kk_flow.angle[i];
                flow.point[i][0] = kk_flow.point[i][0];
                flow.point[i][1] = kk_flow.point[i][1];
        }
        flow.tot_points = kk_flow.tot_points;


        float *FOE = least_square(flow);
        error_distances = dist2lines(flow,FOE[0],FOE[1]);
        //display_flow(flow,FOE);
        //printf("\n dist : %f ",error_distances[flow.tot_points-1]);
        //printf("\n Tot lines = %d",flow.tot_points);
        //printf(", mean error : %f ",mean_err);
        //printf(", FOE = %f, %f",FOE[0], FOE[1]);

        kk = kk + 1;

    }

    // Second implementation of RANSAC
    float  best_inters = INFINITY;
    opt_flow best_flow;
    for(int i=0;i<500;i++)
    {
        opt_flow inliers = random_sample(flow,200);
        float *FOE = least_square(inliers);
        float *error_dist = dist2lines(inliers,FOE[0],FOE[1]);
        float error_inters = 0;
        for(int j=0;j<inliers.tot_points;j++)
            error_inters = error_inters + error_dist[j];
        if(error_inters<best_inters)
        {
            /*
                Copy inliers to best_flow
            */
            best_flow.flow_mag = (float *)malloc(sizeof(float));
            best_flow.angle =    (float *)malloc(sizeof(float));
            best_flow.tot_points = 0;
            for(int i=0;i<inliers.tot_points;i++)
            {
                // Make space and copy the object
                best_flow.flow_mag = (float *)realloc(best_flow.flow_mag,(best_flow.tot_points+1)*sizeof(float));
                best_flow.angle = (float *)realloc(best_flow.angle,(best_flow.tot_points+1)*sizeof(float));
                best_flow.point = (int **)realloc(best_flow.point,(best_flow.tot_points+1)*sizeof(*best_flow.point));
                best_flow.point[best_flow.tot_points] = (int*)malloc(2*sizeof(int));

                best_flow.flow_mag[best_flow.tot_points] = inliers.flow_mag[i];
                best_flow.angle[best_flow.tot_points] = inliers.angle[i];
                best_flow.point[best_flow.tot_points][0] = inliers.point[i][0];
                best_flow.point[best_flow.tot_points][1] = inliers.point[i][1];
                best_flow.tot_points = best_flow.tot_points + 1;
            }

            best_inters = error_inters;
            //printf("\n best_inters = %f ",best_inters);
        }
    }

    return best_flow;
}

/*
    Function to randomly sample n flow vectors
*/
opt_flow random_sample(opt_flow flow, int n)
{
    if(n>flow.tot_points)           // in case of limited flow vectors
    {
        printf("\n Full sampling inside RANSAC! ");
        n = flow.tot_points;
    }

    int arr[flow.tot_points];
    for (int i = 0; i < flow.tot_points; i++)
            arr[i] = i;
    // shuffle array
    for (int i = 0; i < flow.tot_points; i++)
        {
            int temp = arr[i];
            int randomIndex = rand() % flow.tot_points;

            arr[i] = arr[randomIndex];
            arr[randomIndex] = temp;
        }
    // select first n sample from arr
    int sample[n];
    //printf("\n");
    for (int i = 0; i < n; i++)
        {
            sample[i] = arr[i];
            //printf(" %d ",sample[i]);
        }

    opt_flow flow_sample;
    flow_sample.flow_mag = (float *)malloc(sizeof(float));
    flow_sample.angle =    (float *)malloc(sizeof(float));
    flow_sample.tot_points = 0;
    for(int i=0;i<n;i++)
    {
        flow_sample.flow_mag = (float *)realloc(flow_sample.flow_mag,(flow_sample.tot_points+1)*sizeof(float));
        flow_sample.angle = (float *)realloc(flow_sample.angle,(flow_sample.tot_points+1)*sizeof(float));
        flow_sample.point = (int **)realloc(flow_sample.point,(flow_sample.tot_points+1)*sizeof(*flow_sample.point));
        flow_sample.point[flow_sample.tot_points] = (int*)malloc(2*sizeof(int));

        flow_sample.flow_mag[flow_sample.tot_points] = flow.flow_mag[sample[i]];
        flow_sample.angle[flow_sample.tot_points] = flow.angle[sample[i]];
        flow_sample.point[flow_sample.tot_points][0] = flow.point[sample[i]][0];
        flow_sample.point[flow_sample.tot_points][1] = flow.point[sample[i]][1];
        flow_sample.tot_points = flow_sample.tot_points +1;
    }

    return flow_sample;
}

/*
    Remove the flow vectors that falls outside the threshold bounds for the standard deviation
    of the error distances,ie. remove the outlier lines
*/
opt_flow remove_outlier_flow(opt_flow flow, float* error_distances,float mean_err,float err_std_thresh)
{
        float std_err = std_dev(error_distances,flow.tot_points);
        static opt_flow new_flow;
        new_flow.flow_mag = (float *)malloc(sizeof(float));
        new_flow.angle =    (float *)malloc(sizeof(float));
        new_flow.tot_points = 0;
        for(int i=0;i<flow.tot_points;i++)
        {
            if(error_distances[i]>=mean_err-err_std_thresh*std_err && \
               error_distances[i]<=mean_err+err_std_thresh*std_err)
            {
                // Make space and copy the object
                new_flow.flow_mag = (float *)realloc(new_flow.flow_mag,(new_flow.tot_points+1)*sizeof(float));
                new_flow.angle = (float *)realloc(new_flow.angle,(new_flow.tot_points+1)*sizeof(float));
                new_flow.point = (int **)realloc(new_flow.point,(new_flow.tot_points+1)*sizeof(*new_flow.point));
                new_flow.point[new_flow.tot_points] = (int*)malloc(2*sizeof(int));

                new_flow.flow_mag[new_flow.tot_points] = flow.flow_mag[i];
                new_flow.angle[new_flow.tot_points] = flow.angle[i];
                new_flow.point[new_flow.tot_points][0] = flow.point[i][0];
                new_flow.point[new_flow.tot_points][1] = flow.point[i][1];
                new_flow.tot_points = new_flow.tot_points + 1;
            }
        }
        //printf("\n tot_points : %d ",new_flow.tot_points);
        return new_flow;
}

/*
    Evaluate mean of n data elements
*/
float mean(float* a, int n)
{
    float sum=0;
    for(int i=0;i<n;i++)
        sum = sum + a[i];
    return sum/(float)n;
}

/*
    Evaluate standard deviation of n data elements
*/
float std_dev(float* a, int n)
{
    float sum = 0, avg = mean(a,n);
    for(int i=0;i<n;i++)
        sum = sum + pow(a[i]-avg,2);
    return sqrt(sum/(float)n);
}

/*
    Evaluate perpendicular distances from FOE to all of the lines
*/
float* dist2lines(opt_flow flow, float FOEx, float FOEy)
{
    //printf("\n FOE = %f, %f",FOEx, FOEy);
    float *error_distances = (float *)malloc(flow.tot_points*sizeof(float));
    //float error_distances[flow.tot_points];
    for(int i=0;i<flow.tot_points;i++)
    {
        error_distances[i] = abs((tan(flow.angle[i]))*FOEx-FOEy+flow.point[i][1]- \
                                 (tan(flow.angle[i]))*flow.point[i][0])/sqrt(1+pow((tan(flow.angle[i])),2));
        //printf("\n Errors  = %f",error_distances[i]);
    }
    return error_distances;
}

/*
    Filter out those flow vectors with zero flow magnitude in the direction perpendicular
    to the edges
*/
opt_flow find_no_flow(opt_flow flow)
{
    // keep low flow magnitude edge points and remove singularity points

    float flow_mag_thresh = 1;
    opt_flow no_flow;
    no_flow.flow_mag = (float *)malloc(sizeof(float));
    no_flow.angle =    (float *)malloc(sizeof(float));
    no_flow.tot_points = 0;
    for(int i=0;i<flow.tot_points;i++)
    {
        //printf("\n sdd %f ",(round((flow.angle[i]-PI/2)*10)/10));
        if(flow.flow_mag[i]<flow_mag_thresh && (round((flow.angle[i]-PI/2)*10)/10)!=0.0&& \
                                               (round((flow.angle[i]+PI/2)*10)/10)!=0.0)
        //&& (round((flow.angle[i])*10)/10)!=0.0 && (round((flow.angle[i])*10)/10)!=3.1 && (round((flow.angle[i])*10)/10)!=-3.1)
        {
            // Make space and copy the object
            //printf("\n %f ", (round((flow.angle[i]-PI/2)*10)/10));
            no_flow.flow_mag = (float *)realloc(no_flow.flow_mag,(no_flow.tot_points+1)*sizeof(float));
            no_flow.angle = (float *)realloc(no_flow.angle,(no_flow.tot_points+1)*sizeof(float));
            no_flow.point = (int **)realloc(no_flow.point,(no_flow.tot_points+1)*sizeof(*no_flow.point));
            no_flow.point[no_flow.tot_points] = (int*)malloc(2*sizeof(int));
            //printf("\n mag %f ",flow.flow_mag[i]);

            no_flow.flow_mag[no_flow.tot_points] = flow.flow_mag[i];
            no_flow.angle[no_flow.tot_points] = flow.angle[i];
            no_flow.point[no_flow.tot_points][0] = flow.point[i][0];
            no_flow.point[no_flow.tot_points][1] = flow.point[i][1];
            no_flow.tot_points = no_flow.tot_points +1;
        }
    }
    //printf("\n tot_points %d ",no_flow.tot_points);
//    for (int i=0; i<no_flow.tot_points; i++)
//    {
//        printf("\n angle %f ",no_flow.angle[i]);
//    }
    return no_flow;
}

/*
    Given 2 frames, evaluate the optical flow at all edge points along the direction
    perpendicular to the edges
*/
opt_flow get_optical_flow_edges(IplImage *I1,IplImage *I2,bool graphics)
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

//    int kki = 200;
//    for (int kkj=0;kkj<500;kkj++)
//    {schar dxx = CV_IMAGE_ELEM(dx, schar,kki,2*kkj);
//    schar dyy = CV_IMAGE_ELEM(dy, schar,kki,2*kkj);
//    printf("dx=%d, dy=%d, angle=%f \n",dxx,dyy,Gdir[kki][kkj]);
//    }

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
//                printf("\n\npoint=[%d,%d], ang=%f, xb= ",point[tot_points][0],point[tot_points][1],\
//                       angle[tot_points]*180/PI);
//                for(int i=0;i<range+1;i++)
//                    printf("%d ",xb[i]);
//                printf(", yb= ");
//                for(int i=0;i<range+1;i++)
//                    printf("%d ",yb[i]);
//                printf(", xt= ");
//                for(int i=0;i<range+1;i++)
//                    printf("%d ",xt[i]);
//                printf(", yt= ");
//                for(int i=0;i<range+1;i++)
//                    printf("%d ",yt[i]);
                }

                // SSD Correlation along the line and evaluate optical flow

                int window = 15,best_match[2] = {0,0};
                int best_dir = 1, best_i = 0;
                unsigned long int best_corr = INFINITY, all_corr[2*range+1];
                for(int i=0;i<2*range+1;i++)
                    all_corr[i] = 0;

                for(int i=0;i<range+1;i++)
                {
                    // check if the search boundaries doesn't go out of the image
                    if (yb[i]-(window-1)/2>=0 && yb[i]+(window-1)/2<siz.height && \
                        xb[i]-(window-1)/2>=0 && xb[i]+(window-1)/2<siz.width && \
                        yt[i]-(window-1)/2>=0 && yt[i]+(window-1)/2<siz.height && \
                        xt[i]-(window-1)/2>=0 && xt[i]+(window-1)/2<siz.width && \
                        point[tot_points][0]-(window-1)/2>=0 && point[tot_points][0]+(window-1)/2<siz.height && \
                        point[tot_points][1]-(window-1)/2>=0 && point[tot_points][1]+(window-1)/2<siz.width)
                    {
                        unsigned long int corr_b, corr_t, corr_d;
                        bool if_ssd = true;
                        int corr_match[2], corr_dir;
                        corr_b = cross_corr(I1_gray,I2_gray,j,k,window,xb[i],yb[i],if_ssd);
                        corr_t = cross_corr(I1_gray,I2_gray,j,k,window,xt[i],yt[i],if_ssd);
                        //printf("\n corr_b %lu, corr_t %lu ",corr_b,corr_t);
                        //printf("\n %u ",CV_IMAGE_ELEM(I2_gray,uchar,yb[i]-(window-1)/2,xb[i]-(window-1)/2)- \
                    CV_IMAGE_ELEM(I1_gray,uchar,j-(window-1)/2,k-(window-1)/2));

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
                        //printf("\n %lu ",all_corr[range+i]);

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
                        //printf("\n corr_d %lu ",corr_d);

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

                //printf("\n best_corr %lu ",best_corr);
//                if(best_corr!=4294967295)
//                {
//                    //printf("\n best_corr %li ",best_corr);
//                printf("\n\n All corr = ");
//                for(int ii=0;ii<2*range+1;ii++)
//                {
//                    printf(" %lu ",all_corr[ii]);
//                }
//                }


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
                if (flow_mag[tot_points-1]>flow_thresh || best_corr==4294967295)
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
                //printf("\n i = %d, flow = %f",i+1,flow_mag[i]);
                //cvLine(I1_gray,cvPoint(100,100),cvPoint(200,200),cvScalar(0,0,255),2,8,0);
            }
        //cvLine(I1_gray,cvPoint(100,100),cvPoint(200,200),cvScalar(0,0,255),2,8,0);
        //cvCircle(I1_gray,cvPoint(100,100),5,cvScalar(0,255,0),0);
        printf("\n\n %d",tot_points);
        cvShowImage("new", I1);
        //cvShowImage("image", angle);
        cvWaitKey(0);
    }
//    printf("\n correlation = %lu ",cross_corr(I1_gray,I2_gray,200,200,15, 210,210, true));
//    for(int m=0;m<15;m++)
//    {
//        printf("\n ");
//        for(int n=0;n<15;n++)
//        {
//            printf(" %d ",CV_IMAGE_ELEM(I2_gray,uchar,210-(15-1)/2+m,210-(15-1)/2+n)-\
//                   CV_IMAGE_ELEM(I1_gray,uchar,201-(15-1)/2+m,201-(15-1)/2+n));
//        }
//    }


    opt_flow flow;
    flow.angle = angle;
    flow.flow_mag = flow_mag;
    flow.point = point;
    flow.tot_points = tot_points;
    return flow;
}

/*
   Remove the element at the end of the dynamically allocated array
*/
int** remove_end(int** arr, int sizeOfArray)
{
    // Function to remove the last element of a dynamically allocated 2D array
    int** temp = (int **)malloc((sizeOfArray-1)*sizeof(*arr));
    if (sizeOfArray != 0)
        memcpy(temp, arr, (sizeOfArray - 1) * sizeof(*arr));
    free(arr);
    return temp;
}

/*
   Remove the element at the end of the dynamically allocated array
*/
float* remove_end(float* arr, int sizeOfArray)
{
    // Function to remove the last element of a dynamically allocated array
    float* temp = (float *)malloc((sizeOfArray - 1) * sizeof(float));
    if (sizeOfArray != 0)
        memcpy(temp, arr, (sizeOfArray - 1) * sizeof(float));
    free(arr);
    return temp;
}

/*
   Evaluate cross correlation between 2 image segments
*/
unsigned long int cross_corr(IplImage *I1_gray,IplImage *I2_gray,int j,int k,int window,int x,int y,bool if_ssd)
{
    // Function to evaluate the cross correlation of 2 image segments
    unsigned long int corr = 0;
    if(if_ssd)
    {
        // SSD correlation
        for(int m=0;m<window;m++)
            for(int n=0;n<window;n++)
            {
                schar c_diff = CV_IMAGE_ELEM(I2_gray,uchar,y-(window-1)/2+m,x-(window-1)/2+n)- \
                               CV_IMAGE_ELEM(I1_gray,uchar,j-(window-1)/2+m,k-(window-1)/2+n);
                corr = corr + c_diff*c_diff;
                //printf("\nm = %d, n = %d, c_diff = %d %d, corr = %lu ", m,n,c_diff,CV_IMAGE_ELEM(I2_gray,uchar,y-(window-1)/2+m,x-(window-1)/2+n)- \
                               CV_IMAGE_ELEM(I1_gray,uchar,j-(window-1)/2+m,k-(window-1)/2+n),corr);
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

/*
   Evaluate gradient direction at all edge points
*/
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
            angle[i][j] = atan2(dyy,-dxx)*180/PI;
            //printf("%d  ",dyy);
        }
    }
    return angle;
}

/*
   Find the Sobel edges in the image
*/
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

/*
   Evaluate inverse of a square matrix
*/
float** mat_inv(float **a,int m)
{
    float **C = (float **)malloc(m*sizeof(float *));
    for (int i=0; i<m; i++)
        C[i] = (float *)malloc(m * sizeof(float));
    float **adj = Transpose(CoFactor(a,m),m,m);
    float det = Determinant(a,m);
    for(int i=0;i<m;i++)
        for(int j=0;j<m;j++)
            C[i][j] = adj[i][j]/det;
    return C;
}

/*
   Evaluate Pseudo inverse of a matrix
*/
float** pinv(float **A,int m,int n)
{
    float **At = Transpose(A,m,n);
    float **AtA = mat_multiply(At,n,m,A,m,n);
    return mat_multiply(mat_inv(AtA,n),n,n,At,n,m);
}

/*
   Multiply 2 matrices
*/
float** mat_multiply(float **A,int m, int n, float **B, int o, int p)
{
    if(n!=o)
        printf("\n Matrix cannot be multiplied! ");

    float **C = (float **)malloc(m*sizeof(float *));
    for (int i=0; i<m; i++)
        C[i] = (float *)malloc(p * sizeof(float));

     for(int i=0;i<m;i++) //row of first matrix
        {
            for(int j=0;j<p;j++)
                {  //column of second matrix
                    float sum=0;
                    for(int k=0;k<n;k++)
                    {
                        sum = sum+A[i][k]*B[k][j];
                        C[i][j]=sum;
                    }
                }
        }
     return C;
}

/*
   Recursive definition of determinate using expansion by minors.
*/
float Determinant(float **a,int n)
{
   int i,j,j1,j2;
   float det = 0;
   float **m = NULL;

   if (n < 1) { /* Error */

   } else if (n == 1) { /* Shouldn't get used */
      det = a[0][0];
   } else if (n == 2) {
      det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
   } else {
      det = 0;
      for (j1=0;j1<n;j1++) {
         m = (float **)malloc((n-1)*sizeof(float *));
         for (i=0;i<n-1;i++)
            m[i] = (float *)malloc((n-1)*sizeof(float));
         for (i=1;i<n;i++) {
            j2 = 0;
            for (j=0;j<n;j++) {
               if (j == j1)
                  continue;
               m[i-1][j2] = a[i][j];
               j2++;
            }
         }
         det += pow(-1.0,j1+2.0) * a[0][j1] * Determinant(m,n-1);
         for (i=0;i<n-1;i++)
            free(m[i]);
         free(m);
      }
   }
   return(det);
}

/*
   Find the cofactor matrix of a square matrix
*/
float** CoFactor(float **a,int n)
{
   int i,j,ii,jj,i1,j1;
   float det;
   float **b = (float **)malloc(n*sizeof(float *));
   for (int i=0; i<n; i++)
        b[i] = (float *)malloc(n * sizeof(float));

   float **c;
   c = (float **)malloc((n-1)*sizeof(float *));
   for (i=0;i<n-1;i++)
     c[i] = (float *)malloc((n-1)*sizeof(float));

   for (j=0;j<n;j++) {
      for (i=0;i<n;i++) {

         /* Form the adjoint a_ij */
         i1 = 0;
         for (ii=0;ii<n;ii++) {
            if (ii == i)
               continue;
            j1 = 0;
            for (jj=0;jj<n;jj++) {
               if (jj == j)
                  continue;
               c[i1][j1] = a[ii][jj];
               j1++;
            }
            i1++;
         }

         /* Calculate the determinate */
         det = Determinant(c,n-1);

         /* Fill in the elements of the cofactor */
         b[i][j] = pow(-1.0,i+j+2.0) * det;
      }
   }
   for (i=0;i<n-1;i++)
      free(c[i]);
   free(c);

   return b;
}

/*
   Transpose of a matrix
*/
float** Transpose(float **a,int m, int n)
{
   int i,j;
   float **c = (float **)malloc(n*sizeof(float *));
   for (i=0;i<n;i++)
     c[i] = (float *)malloc(m*sizeof(float));

   for (i=0;i<m;i++)
      for (j=0;j<n;j++)
         c[j][i] = a[i][j];

   return c;
}

/*
    Quiver equivalent in MATLAB
*/
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

/*
    Function to test out stuff
*/
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
