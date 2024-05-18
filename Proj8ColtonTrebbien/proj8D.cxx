
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkDataSetReader.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkLookupTable.h>
#include <vtkContourFilter.h>
#include <vtkDataSetWriter.h>
#include <vtkCutter.h>
#include <vtkPlane.h>
#include <vtkNew.h>
#include <vtkUnsignedCharArray.h>



int main(int argc, char *argv[])
{
   vtkDataSetReader *reader = vtkDataSetReader::New();
   reader->SetFileName("noise.vtk");
   reader->Update();

   /*
   vtkContourFilter *cf = vtkContourFilter::New();
   cf->SetNumberOfContours(2);
   cf->SetValue(0,4);
   cf->SetValue(1,2.4);
   cf->SetInputConnection(reader->GetOutputPort());
   
   vtkDataSetWriter *wrtr = vtkDataSetWriter::New();
   wrtr->SetFileName("contour.vtk");
   wrtr->SetInputConnection(cf->GetOutputPort());
   wrtr->Write();
   */

   vtkPlane *plane = vtkPlane::New();
   plane->SetOrigin(0,0,0);
   plane->SetNormal(0,0,1);

   vtkCutter *cutter = vtkCutter::New();
   cutter->SetCutFunction(plane);
   cutter->SetInputConnection(reader->GetOutputPort());
   cutter->Update();


   vtkDataSetMapper *mapper = vtkDataSetMapper::New();
   mapper->SetInputData(cutter->GetOutput());
   
   vtkLookupTable *lut = vtkLookupTable::New();
   
   //colors->SetNumberOfTuples(256);
   for(unsigned int i = 0; i < 256; i++){
      lut->SetTableValue(i, int(255*i/(256-1)), 0, int(255*(256-i-1)/(256-1)));
   }
   mapper->SetLookupTable(lut);
   mapper->SetScalarRange(1,6);
   lut->Build();
  

   vtkActor *actor = vtkActor::New();
   actor->SetMapper(mapper);

   vtkRenderer *ren = vtkRenderer::New();
   ren->AddActor(actor);

   vtkRenderWindow *renwin = vtkRenderWindow::New();
   renwin->AddRenderer(ren);
   //A)  Make the visualization window be 768x768 when the program is first invoked.   
   renwin->SetSize(768, 768);

   vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
   iren->SetRenderWindow(renwin);
   renwin->Render();
   iren->Start();

  
}


