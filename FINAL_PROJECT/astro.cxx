#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkImageData.h"
#include <vtkPNGWriter.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>
#include <math.h>



/*====================OUTLINE====================

    1. Find Ray For that pixel
    2. interesct volume with ray
    3. calculate color from interestion

        iteration from Fa fback, then those 2 to the next ith iteration fa. 
    4. assign Color to pixel

======================OUTLINE==================*/

struct Camera
{
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];
};


struct TransferFunction
{
    double          min;
    double          max;
    int             numBins;
    unsigned char  *colors;  // size is 3*numBins
    double         *opacities; // size is numBins

    // Take in a value and applies the transfer function.
    // Step #1: figure out which bin "value" lies in.
    // If "min" is 2 and "max" is 4, and there are 10 bins, then
    //   bin 0 = 2->2.2
    //   bin 1 = 2.2->2.4
    //   bin 2 = 2.4->2.6
    //   bin 3 = 2.6->2.8
    //   bin 4 = 2.8->3.0
    //   bin 5 = 3.0->3.2
    //   bin 6 = 3.2->3.4
    //   bin 7 = 3.4->3.6
    //   bin 8 = 3.6->3.8
    //   bin 9 = 3.8->4.0
    // and, for example, a "value" of 3.15 would return the color in bin 5
    // and the opacity at "opacities[5]".
    void ApplyTransferFunction(double value, unsigned char *RGB, double &opacity, int i, int j)
    {
        //cout << "Applying transfer function to " << value << endl;
        if((value > max) or (value < min)){
            if((i == 50) and (j ==50)){
                //cout << "Out of range: no valid bin" << endl;
            }
            RGB[0] = 0;
            RGB[1] = 0;
            RGB[2] = 0;
            return;
        }
   
        double step = (max-min)/numBins;
        
        //double binval =  (value - min)/(max-min);
        //cout << "binval: " << binval  << "step: " << step << endl;
        int bin = 0;
        for(int i=0; i < numBins; i++){
            if (value >= (min + (i*step))){
                bin = i;

            }
        }

        if((i == 50) and (j ==50)){
            //cout << "Mapped to bin: " << bin << endl;
        }
        RGB[0] = colors[3*bin+0];
        RGB[1] = colors[3*bin+1];
        RGB[2] = colors[3*bin+2];
        opacity = opacities[bin];
    }
};

TransferFunction
SetupTransferFunction(void)
{
    int  i;

    TransferFunction rv;
    rv.min = 10;
    rv.max = 15;
    rv.numBins = 256;
    rv.colors = new unsigned char[3*256];
    rv.opacities = new double[256];
    unsigned char charOpacity[256] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13, 14, 14, 14, 14, 14, 14, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 5, 4, 3, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 17, 17, 17, 17, 17, 17, 16, 16, 15, 14, 13, 12, 11, 9, 8, 7, 6, 5, 5, 4, 3, 3, 3, 4, 5, 6, 7, 8, 9, 11, 12, 14, 16, 18, 20, 22, 24, 27, 29, 32, 35, 38, 41, 44, 47, 50, 52, 55, 58, 60, 62, 64, 66, 67, 68, 69, 70, 70, 70, 69, 68, 67, 66, 64, 62, 60, 58, 55, 52, 50, 47, 44, 41, 38, 35, 32, 29, 27, 24, 22, 20, 20, 23, 28, 33, 38, 45, 51, 59, 67, 76, 85, 95, 105, 116, 127, 138, 149, 160, 170, 180, 189, 198, 205, 212, 217, 221, 223, 224, 224, 222, 219, 214, 208, 201, 193, 184, 174, 164, 153, 142, 131, 120, 109, 99, 89, 79, 70, 62, 54, 47, 40, 35, 30, 25, 21, 17, 14, 12, 10, 8, 6, 5, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        };

    for (i = 0 ; i < 256 ; i++)
        rv.opacities[i] = charOpacity[i]/255.0;
    const int numControlPoints = 8;
    unsigned char controlPointColors[numControlPoints*3] = { 
           71, 71, 219, 0, 0, 91, 0, 255, 255, 0, 127, 0, 
           255, 255, 0, 255, 96, 0, 107, 0, 0, 224, 76, 76 
       };
    double controlPointPositions[numControlPoints] = { 0, 0.143, 0.285, 0.429, 0.571, 0.714, 0.857, 1.0 };
    for (i = 0 ; i < numControlPoints-1 ; i++)
    {
        int start = controlPointPositions[i]*rv.numBins;
        int end   = controlPointPositions[i+1]*rv.numBins+1;
        //cerr << "Working on " << i << "/" << i+1 << ", with range " << start << "/" << end << endl;
        if (end >= rv.numBins)
            end = rv.numBins-1;
        for (int j = start ; j <= end ; j++)
        {
            double proportion = (j/(rv.numBins-1.0)-controlPointPositions[i])/(controlPointPositions[i+1]-controlPointPositions[i]);
            if (proportion < 0 || proportion > 1.)
                continue;
            for (int k = 0 ; k < 3 ; k++)
                rv.colors[3*j+k] = proportion*(controlPointColors[3*(i+1)+k]-controlPointColors[3*i+k])
                                 + controlPointColors[3*i+k];
        }
    }    

    return rv;
}


Camera
SetupCamera(void)
{
    Camera rv;
    rv.focus[0] = 0;
    rv.focus[1] = 0;
    rv.focus[2] = 0;
    rv.up[0] = 0;
    rv.up[1] = -1;
    rv.up[2] = 0;
    rv.angle = 30;
    rv.near = 7.5e+7;
    rv.far = 1.4e+8;
    rv.position[0] = -8.25e+7;
    rv.position[1] = -3.45e+7;
    rv.position[2] = 3.35e+7;

    return rv;
}

int GetNumberOfPoints(const int *dims)
{
    // 3D
    return dims[0]*dims[1]*dims[2];
    // 2D
    //return dims[0]*dims[1];
}

// ****************************************************************************
//  Function: GetNumberOfCells
//
//  Arguments:
//
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the number of cells in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfCells(const int *dims)
{
    // 3D
    return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    //return (dims[0]-1)*(dims[1]-1);
}


// ****************************************************************************
//  Function: GetPointIndex
//
//  Arguments:
//      idx:  the logical index of a point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always zero if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z
//            2D data sets would have Z=1
//
//  Returns:  the point index
//
// ****************************************************************************

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    //return idx[1]*dims[0]+idx[0];
}


// ****************************************************************************
//  Function: GetCellIndex
//
//  Arguments:
//      idx:  the logical index of a cell.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always zero if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z
//            2D data sets would have Z=1
//
//  Returns:  the cell index
//
// ****************************************************************************

int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)*idx[0];
    // 2D
    //return idx[1]*(dims[0]-1)+idx[0];
}

// ****************************************************************************
//  Function: GetLogicalPointIndex
//
//  Arguments:
//      idx (output):  the logical index of the point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1] 
//              2 <= idx[2] < dims[2] (or always zero if 2D)
//      pointId:  a number between 0 and (GetNumberOfPoints(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
    idx[0] = pointId%dims[0];
    idx[1] = (pointId/dims[0])%dims[1];
    idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    //idx[0] = pointId%dims[0];
    //idx[1] = pointId/dims[0];
}


// ****************************************************************************
//  Function: GetLogicalCellIndex
//
//  Arguments:
//      idx (output):  the logical index of the cell index.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always zero if 2D)
//      cellId:  a number between 0 and (GetNumberOfCells(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
    // 3D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    //idx[0] = cellId%(dims[0]-1);
    //dx[1] = cellId/(dims[0]-1);
}

double Interpulate(double point_ind, double left, double right, double F1_point, double F2_point){

    double t = (point_ind-left)/(right-left);
    double end = F1_point + (t*(F2_point - F1_point));    
    return end;
}


double EvaluateFieldAtLocation(const double *pt, const int *dims, const float *X, 
                              const float *Y, const float *Z, const float *F)
{
    
     if(( pt[0] > X[dims[0]-1]) or( pt[0] < X[0]) or ( pt[1] < Y[0]) or (pt[1] > Y[dims[1]-1]) or( pt[2] > Z[dims[2]-1]) or(pt[2] < Z[0])) {
        return 0;
    }


    double ptx = pt[0];
    double pty = pt[1];  
    double ptz = pt[2];  
    //int cell = GetCellIndex(pt, dims)
    //cout << cell << end;
    double close_x = 0;
    double close_y = 0;
    double close_z = 0;
    int count_x = -1;
    int count_y = -1;
    int count_z = -1;


    //First Find the cell its in.

    for( int i = 0; i < dims[0]; i++){
        if (X[i] < pt[0]){
            close_x = X[i];
            count_x+=1;
        }
        else{
            break;
        }
    }
     for( int i = 0; i < dims[1]; i++){
        if (Y[i] < pt[1]){
            close_y = Y[i];
            count_y+=1; 
        }
        else{
            break;
        }
    }
    for( int i = 0; i < dims[2]; i++){
        if (Y[i] < pt[2]){
            close_z = Z[i];
            count_z+=1; 
        }
        else{
            break;
        }
    }

    int offsetsI[8] = { 0, 1, 1, 0, 0, 1, 1, 0 };
    int offsetsJ[8] = { 0, 0, 0, 0, 1, 1, 1, 1 };
    int offsetsK[8] = { 0, 0, 1, 1, 0, 0, 1, 1 };



    int ind[3];
    ind[0] = count_x;
    ind[1] = count_y;
    ind[2] = count_z;
    //cout << ind[0] << " " << ind[1] << " " << ind[2] << " " << endl;
    //GetLogicalCellIndex(ind, i, dims); 
    int points[8][3];
    //Very Useless Table that should have been done in a loop lmao
    points[0][0] = ind[0] + offsetsI[0];   points[0][1] = ind[1] + offsetsJ[0];   points[0][2] = ind[2] + offsetsK[0]; //V0 Bottom Left front
    points[1][0] = ind[0] + offsetsI[1];   points[1][1] = ind[1] + offsetsJ[1];   points[1][2] = ind[2] + offsetsK[1]; //V1 bottom right front
    points[2][0] = ind[0] + offsetsI[2];   points[2][1] = ind[1] + offsetsJ[2];   points[2][2] = ind[2] + offsetsK[2]; //V2 bottom right back
    points[3][0] = ind[0] + offsetsI[3];   points[3][1] = ind[1] + offsetsJ[3];   points[3][2] = ind[2] + offsetsK[3]; //V3 bottom left back
    points[4][0] = ind[0] + offsetsI[4];   points[4][1] = ind[1] + offsetsJ[4];   points[4][2] = ind[2] + offsetsK[4]; //V4 top left front
    points[5][0] = ind[0] + offsetsI[5];   points[5][1] = ind[1] + offsetsJ[5];   points[5][2] = ind[2] + offsetsK[5]; //v5 top right fron
    points[6][0] = ind[0] + offsetsI[6];   points[6][1] = ind[1] + offsetsJ[6];   points[6][2] = ind[2] + offsetsK[6]; //V6 top right back
    points[7][0] = ind[0] + offsetsI[7];   points[7][1] = ind[1] + offsetsJ[7];   points[7][2] = ind[2] + offsetsK[7]; //V7 top back left

    int point_index_[8];
    //float logical_index[8];
    double verts[8];
 
    for(int x =0; x <8; x++){
        point_index_[x] = GetPointIndex(points[x], dims);
        //cout<< x <<":" << point_index_[x]<<endl;

        verts[x]  = F[point_index_[x]];
        //cout << "vert: " << verts[x] << endl;
    }

     //cout << endl;

    //cout << "points" << points[0][0] << " "<<  points[0][1] << " " << points[0][2] << endl;
    //cout << "find" << count_x << " "<<  count_y << " " << count_z << endl;
    
    double front_bot = Interpulate(pt[0], X[points[0][0]], X[points[1][0]], verts[0], verts[1]);
    double front_top = Interpulate(pt[0], X[points[4][0]], X[points[5][0]], verts[4], verts[5]);
    //cout << "Front Top/Fbot: " << front_top << " " << front_bot <<endl;
    double back_bot  = Interpulate(pt[0], X[points[3][0]], X[points[2][0]], verts[3], verts[2]);
    double back_top  = Interpulate(pt[0], X[points[7][0]], X[points[6][0]], verts[7], verts[6]);
    //cout << "Back Top/Fbot: " << back_top << " " << back_bot <<endl;
    double front = Interpulate(pt[1], Y[points[0][1]], Y[points[4][1]], front_bot, front_top);
    //cout << "front: "<< front << endl;
    double back = Interpulate(pt[1], Y[points[0][1]], Y[points[4][1]], back_bot, back_top);
    // cout << "back: "<< back << endl;
    double final = Interpulate(pt[2], Z[points[0][2]], Z[points[6][2]],front ,back);
    //double font_top =  nterpulate(pt[0], F[points[0][0], F[]])
    //cout << "old_front: " << front_bot << endl;
    // 
    
    return final;
}


void crossProduct(double vect_A[], double vect_B[], double cross_P[])
    //I found this off geeks for geeks 
{
    cross_P[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1];
    cross_P[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2];
    cross_P[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0];
}


void Composite_sample(double *RGB_front, double *RGB_back){

    RGB_back[0] = (RGB_back[0] + (1-RGB_back[3]) * RGB_front[0]/255 * RGB_front[3]); 
    RGB_back[1] = (RGB_back[1] + (1-RGB_back[3]) * RGB_front[1]/255 * RGB_front[3]); 
    RGB_back[2] = (RGB_back[2] + (1-RGB_back[3]) * RGB_front[2]/255 * RGB_front[3]); 
    RGB_back[3] =  RGB_front[3] + ((1-RGB_front[3]) * RGB_back[3]);
}

void
WriteImage(vtkImageData *img, const char *filename)
{
    std::string full_filename = filename;
    full_filename += ".png";
    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(full_filename.c_str());
    writer->Write();
    writer->Delete();
}

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *image = vtkImageData::New();
    image->SetDimensions(width, height, 1);
    //image->SetWholeExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetUpdateExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetNumberOfScalarComponents(3);
    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
    //image->AllocateScalars();

    return image;
}

void ApplyColorMap(double *XRGB, unsigned char *RGB)
{
 
    RGB[0] = int(XRGB[0]*255);

    RGB[1] = int(XRGB[1]*255);
    //RGB[2] = 255 * scale;
    RGB[2] = int(XRGB[2]*255); 
  
}

int main()
{
    TransferFunction tf = SetupTransferFunction();
    
    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("astro512.vtk");
    rdr->Update();
    if (rdr->GetOutput() == NULL || rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Could not find input file." << endl;
        exit(EXIT_FAILURE);
    }
   // cout << "2" << endl;
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    int dims[3];
   // cout << "3" << endl;
    rgrid->GetDimensions(dims);
   // cout << "4" << endl;
    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    //cout << "5" << endl;
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    //cout << "6" << endl;
    float *Z = (float *) rgrid->GetZCoordinates()->GetVoidPointer(0);
    //cout << "7" << endl;
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    //cerr << "8" << endl;
    int pixels = 512;
    int SPR = 1000;

    int nx = 512;
    int ny = 512;
    //cerr << "made it" << endl;

    vtkImageData *image;
    unsigned char *buffer;

    image = NewImage(nx, ny);
    buffer = (unsigned char *) image->GetScalarPointer(0,0,0);


    for (int j = 0 ; j < 3*nx*ny ; j++){
            buffer[j] = 0;
    }

    char file[13] = "astro512.vtk";
    cout << "Program specs:" << endl;
    cout << endl;
    cout << "IMAGE SIZE: " << pixels << "x" << pixels << endl;
    cout << endl;
    cout << "SAMPLES PER RAY: " << SPR << endl;
    cout << endl;
    cout << "Data input: " << file << endl; 
    cout << endl;
    cout << "All print statements are for the ray corresponding to pixel (50,50)" << endl;
    cout << endl;
    cout << endl;

    Camera my_camera;
    my_camera = SetupCamera();
    //cout << my_camera.near << endl;

    //MATH PART TO CODE... 
    //Start with the focol point. (that should be 50x50 pixel)
    //Focus - position gives you the LOOK vector.
    //fov x and y are just angle. its just 1 angle, the value of angle g.
    //step size, 

    double look[3];
    look[0] = my_camera.focus[0] - my_camera.position[0];
    look[1] = my_camera.focus[1] - my_camera.position[1];
    look[2] = my_camera.focus[2] - my_camera.position[2];
   
    //Set Up Our vectors
    double u[3];
    double v[3];
    double top[3];
    double bot[3];
    double delta_x[3];
    double delta_y[3];
    double RayDir[3];
    //double delta_y[3];
    /*Get U vector, by crossproduct look x up
                                    ---------
                                   |look x up|  

    and divison by the magnitutde of the crossproduct sqrt(pow(2)) */
    crossProduct(look, my_camera.up,top);
    double mag = sqrt(pow(top[0],2) + pow(top[1],2) + pow(top[2],2));
 
    u[0] = top[0]/mag;
    u[1] = top[1]/mag;
    u[2] = top[2]/mag;
    /*Get V vector, by crossproduct look x up
                                    ---------
                                   |look x up|  

    and divison by the magnitutde of the crossproduct sqrt(pow(2)) */
    crossProduct(look, u, bot);
    double mag_2 = sqrt(pow(bot[0],2) + pow(bot[1],2) + pow(bot[2],2));
    v[0] = bot[0]/mag_2;
    v[1] = bot[1]/mag_2;
    v[2] = bot[2]/mag_2;

    double tangent = tan(my_camera.angle/2*3.1415926/180.0);
    double math = (2.0 * tangent)/ pixels;

    delta_x[0] = math * u[0];
    delta_x[1] = math * u[1];
    delta_x[2] = math * u[2];

    delta_y[0] = math * v[0];
    delta_y[1] = math * v[1];
    delta_y[2] = math * v[2];

    cout << "R_u: " << u[0] << "," << u[1] << "," << u[2] << endl;
    cout << "R_v: " << v[0] << "," << v[1] << "," << v[2] << endl;
    cout << "R_x: " << delta_x[0] << ", " << delta_x[1] << ", " << delta_x[2] << endl;
    cout << "R_y: " << delta_y[0] << ", " << delta_y[1] << ", " << delta_y[2] << endl;
    cout << "Ray origin = "<< look[0] << ", " << look[1] << ", " << look[2] << endl;

    double lookmag = sqrt(pow(look[0],2) + pow(look[1],2) + pow(look[2],2));

    double a = 0;
    double b = 0;

    double samples = 1000;
    double sample_pos[3];

    double dist = 0.0;

    double stepsize = (my_camera.far - my_camera.near)/samples;
    
    double value = 0.0;
    double values[int(samples)];

    unsigned char RGB[4];
    double DRGB[4];
    double RGB_run[4] = {0.0,0.0,0.0,0.0};
    
    double opacity=0;
    double A=0;
    //cerr << "here" << endl;
   
    for(int i = 0; i < pixels; i++){
        //cerr << "pixel: " << i <<endl;
        for(int j = 0; j < pixels; j++){
            //cerr << "pixel: " << i << ", " << j << endl;

            int offset = 3*(j*nx+i);

            for(int x = 0; x < 3; x++){
                RayDir[x] = (look[x]/lookmag) + ((((2*i)+1-pixels)/2.0) * delta_x[x]) + (((2*j+1-pixels)/2.0) *delta_y[x]); 
            } 
            //cerr << "here" << endl;
            if( (i == 50) and (j ==50)){
                //cout << "Ray direction = "<< RayDir[0] << ", " << RayDir[1] << ", " << RayDir[2] << endl;
                //cout << "Step size is "  << stepsize << endl;
            }

            for(int sam = 0; sam < samples; sam++){
                dist = my_camera.near + (stepsize*sam) ;
                //cerr << "SAM:" << sam << endl;
                for(int x = 0; x < 3; x++){   
                    sample_pos[x] = my_camera.position[x] + (dist * RayDir[x]);
                }
                value = EvaluateFieldAtLocation(sample_pos, dims, X, Y, Z, F);
                values[sam] = value;

                 //cerr << "here" << endl;
                if((i == 50) and (j ==50)){
                    //cout << "Position for sample["<< sam << "] = " << sample_pos[0] << ", " <<  sample_pos[1]  << ", " << sample_pos[2]  << endl;
                    //cout << "value at that postion is " << values[sam] << endl;
                }
            }
            for(int sam = 0; sam < samples; sam++){
                dist = my_camera.near + (stepsize*sam) ;
    
                for(int x = 0; x < 3; x++){   
                    sample_pos[x] = my_camera.position[x] + (dist * RayDir[x]);
                }

                tf.ApplyTransferFunction(values[sam], RGB, opacity,i,j);
                A = 1 - pow((1 - opacity), (500/samples));
                DRGB[0] = (int) RGB[0];
                DRGB[1] = (int) RGB[1];
                DRGB[2] = (int) RGB[2];
                DRGB[3] = A;
                Composite_sample(DRGB, RGB_run);

                if((i == 50) and (j ==50)){
                    //cout << "value at that postion is " << value << endl;
                    //cout << "Applying transfer function to " << values[sam] << endl;
                    //tf.ApplyTransferFunction(values[sam], RGB, opacity,i,j);
                   // cout << "Sample[" << sam << "] was mapped by the transfer function to " << (int)RGB[0] << "," << (int)RGB[1] << "," << (int)RGB[2] <<", opacity= " << opacity << endl;
                    
                   // double A = 1 - pow((1 - opacity), (500/samples));
                    //cout << "After opacity correction, opacity is " << A << endl;
                    /*
                    DRGB[0] = (int) RGB[0];
                    DRGB[1] = (int) RGB[1];
                    DRGB[2] = (int) RGB[2];
                    DRGB[3] = A;
                    */
                    //cout << "R:"<<  DRGB[0] << ", " << "G:"<<  DRGB[1] << ", " <<  "B:" <<  DRGB[2] << ", "<< "A:"<<  DRGB[3] <<  endl;

                    //Composite_sample(DRGB, RGB_run);
                    //cout << "After compositing sample " << sam<< ", the running values are R:"<< RGB_run[0]<< ", " << "G:"<< RGB_run[1]<< ", " <<  "B:" << RGB_run[2] << ", "<< "A:"<< RGB_run[3]<<  endl;   
                }
                
                //cout << "Final color for pixel is " << int(RGB_run[0] *255) << " " << int(RGB_run[1] *255) << " " << int(RGB_run[2] *255) << endl;
                
                //ApplyColorMap(RGB_run, buffer+offset);
            }
        ApplyColorMap(RGB_run, buffer+offset);
            
        //simple fix to reset numbers
        for(int b = 0; b < 4; b++){
                RGB_run[b] = 0;
        }
        A=0;
        opacity= 0;    
        }
        //cerr << "pixel: " << i << endl;
        //cerr << "pixel: " << i << ", " << j << endl;
    }

    
    WriteImage(image, "astro512");
 
}
 
