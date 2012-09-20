#include <iostream>
#include <fstream>
#include <limits>
#include <vector>
#include <sstream>

#define MINI(X,Y) ( (X) < (Y) ? (X) : (Y) )
#define MAXI(X,Y) ( (X) > (Y) ? (X) : (Y) )

#include "CImg.h"

using namespace cimg_library;

/*

 Compile with :
 g++ *.cpp -o hyplot -lm -lpthread -lX11

*/




int main(int argc, char** argv)
{
  int           dimX           = 500;
  int           dimY           = 500;
  int           dimLegendX     = 100;
  int           margeX         = 0;
  int           margeY         = 0;
  int           colorBG        = 0;
  std::string   title          = "Image";
  std::string   imgName        = "Image";
  bool          drawAxe        = false;
  std::string   legend         = "";
  
  unsigned int  nLine          = 1;
  
  bool          isAxeSet       = false;
  bool          isColorBGSet   = false;
  bool          isDataTypeSet  = false;
  bool          dataTypeXYZ    = false;
  bool          isDimSet       = false;
  bool          isTitleSet     = false;
  bool          isMinMaxXSet   = false;
  bool          isMinMaxYSet   = false;
  bool          isMinMaxZSet   = false;
  bool          isLegendSet    = false;
  bool          isDisplaySet   = false;
  bool          display        = false;

  double minX =  std::numeric_limits<double>::max();
  double maxX = -std::numeric_limits<double>::max();
  double minY =  std::numeric_limits<double>::max();
  double maxY = -std::numeric_limits<double>::max();
  double minZ =  std::numeric_limits<double>::max();
  double maxZ = -std::numeric_limits<double>::max();

  if(argc > 1)
  {
    // input file
    std::string fileName = argv[1];
    
    // output picture (bmp file)
    std::string outImgName = argv[1];
    outImgName.append(".jpg");
    
    // output log file
    std::string fileLogName = std::string(argv[1]).append(".log");

    // Opening files
    std::fstream fileLog(fileLogName.c_str(),std::ios::out);
    std::fstream file(argv[1],std::ios::in);
    if(file)
    {
      std::string line;
      std::string keyword;
      std::stringstream sline;
      bool read_header = true;
      double x, y, z;   // coordonate
      int r=0;          // color
      int g=0;
      int b=0;
      
      ////////////////////////////////////////////////////////////////////
      // Reads the header
      ////////////////////////////////////////////////////////////////////
      
      while(read_header && getline(file,line))
      {
        if(line != "") // avoids empty lines
        {
          sline.clear();
          sline.str("");
          sline << line;
          sline >> keyword;

          // Parser
          if(keyword[0] == '#')
          {
            if      (keyword == "#AXE" && !isAxeSet)
            {
              sline >> keyword;
              drawAxe = (keyword == "YES");
              isAxeSet = true;
            }
            else if (keyword == "#BEGIN")
            {
              read_header = false;
            }
            else if (keyword == "#COLOR_BG" && !isColorBGSet)
            {
              sline >> colorBG;
              isColorBGSet = true;
            }
            else if (keyword == "#DATA" && !isDataTypeSet)
            {
              sline >> keyword;
              dataTypeXYZ   = (keyword == "XYZ");
              isDataTypeSet = true;
            }
            else if (keyword == "#DIM" && !isDimSet)
            {
              int dimXTmp, dimYTmp;
              
              sline >> dimXTmp >> dimYTmp;
              if(dimXTmp < 20 || dimYTmp < 20)
              {
                std::cerr << "Warning : the given dimensions are too small. Set to (500x500)." << std::endl;
              }
              else
              {
                dimX = dimXTmp;
                dimY = dimYTmp;
              }
              isDimSet = true;
            }
            else if (keyword == "#DISPLAY" && !isDisplaySet)
            {
              sline >> keyword;
              display = (keyword == "YES");
            }
            else if (keyword == "#LEGEND" && !isLegendSet)
            {
              sline >> legend;
              isDataTypeSet = true;
            }
            else if(keyword == "#MINMAXX" && !isMinMaxXSet)
            {
              sline >> minX >> maxX;
              isMinMaxXSet = true;
            }
            else if(keyword == "#MINMAXY" && !isMinMaxYSet)
            {
              sline >> minY >> maxY;
              isMinMaxYSet = true;
            }
            else if(keyword == "#MINMAXZ" && !isMinMaxZSet)
            {
              sline >> minZ >> maxZ;
              isMinMaxZSet = true;
            }
            else if (keyword == "#TITLE" && !isTitleSet)
            {
              sline >> title;
              imgName = title;
              title =  "Hyplot-   "+title;
              isTitleSet = true;
            }
            
          }
          else
          {
            //
          }
        } // (line != "")
        
      ++nLine;
      }

      ////////////////////////////////////////////////////////////////////
      // Recaps the option
      ////////////////////////////////////////////////////////////////////

      fileLog << " ---------------------------------- " << std::endl;
      fileLog << "---------     HYPLOT-     ----------" << std::endl;
      fileLog << "------------------------------------" << std::endl << std::endl;
      fileLog << "Title                   : " << title << std::endl;
      fileLog << "Dimensions              : " << dimX << "x" << dimY << std::endl;
      fileLog << "Data type               : ";
      if(dataTypeXYZ)
        fileLog << "(X,Y,Z)";
      else
        fileLog << "(X,Y)";
      fileLog << std::endl;
      fileLog << "Color of the background : (" << colorBG << ")" << std::endl;
      fileLog << "Axes will ";
      if(!drawAxe)
        fileLog << "not ";
      fileLog << "be drawn" << std::endl;
      fileLog << " ---------------------------------- " << std::endl;

      ////////////////////////////////////////////////////////////////////
      // Prepares the image with the given options
      ////////////////////////////////////////////////////////////////////
      
      CImg<> img(dimX + 2*margeX, dimY + 2*margeY,1,3,0);
      img.fill(colorBG);
      std::fstream file(argv[1],std::ios::in);
      

      ////////////////////////////////////////////////////////////////////
      // Stores the data, counts the number of lines, searches the minima
      ////////////////////////////////////////////////////////////////////
     
      std::vector<double> vectX;
      std::vector<double> vectY;
      std::vector<double> vectZ;
      unsigned int card=0;       // number of data
      
      std::cout << "Preparing the data..." << std::endl;
      while(getline(file,line))
      {
        sline.clear();
        sline.str("");
        sline << line;
        sline >> x >> y;

        if(dataTypeXYZ)
        {
          sline >> z;
        }
        vectX.push_back(x);
        vectY.push_back(y);
        vectZ.push_back(z);
        
        if(!isMinMaxXSet)
        {
          minX = MINI (minX,x);
          maxX = MAXI (maxX,x);
        }
        if(!isMinMaxYSet)
        {
          minY = MINI (minY,y);
          maxY = MAXI (maxY,y);
        }
        if(!isMinMaxZSet)
        {
          minZ = MINI (minZ,z);
          maxZ = MAXI (maxZ,z);
        }
        
        ++card;
      }

      fileLog << "Number of lines : " << card << std::endl;
      fileLog << "minX : " << minX << std::endl;
      fileLog << "maxX : " << maxX << std::endl;
      fileLog << "minY : " << minY << std::endl;
      fileLog << "maxY : " << maxY << std::endl;
      fileLog << "minZ : " << minZ << std::endl;
      fileLog << "maxZ : " << maxZ << std::endl;

      
      ////////////////////////////////////////////////////////////////////
      // Plots the data
      ////////////////////////////////////////////////////////////////////
      std::cout << "Ploting..." << std::endl;
      
      double coefX;
      if(maxX != minX)
        coefX = ((double)(dimX))/(maxX-minX);
      else
        coefX = 0;
      double offsetX = -coefX*minX;
      
      double coefY;
      if(maxY != minY)
        coefY = ((double)(dimY))/(maxY-minY);
      else coefY = 0;
      double offsetY = -coefY*minY;
      
      double coefZ;
      if(maxZ != minZ)
        coefZ = (255.0)/(maxZ-minZ);
      else coefZ = 0;
      double offsetZ = -coefZ*minZ;
      
      for(int i=0;i<vectX.size();++i)
      {
        printf("%.2lf%%\r",(double)(i*100)/vectX.size());
        b = (int)(vectZ[i]*coefZ + offsetZ);
        r = 255 - b;
        g = r;
        int bb = r;
        const unsigned char color[] = {r,g,bb};
        img.draw_point(margeX + vectX[i]*coefX + offsetX, dimY + margeY - vectY[i]*coefY - offsetY,color);
      }
      
      ////////////////////////////////////////////////////////////////////
      // Legend
      ////////////////////////////////////////////////////////////////////
      
      const unsigned char colorBackGround[] = {colorBG,0,0};
      const unsigned char colorLegend[]     = {255-colorBG,255,255};
      img.draw_text(dimX + 3*margeX, margeY, legend.c_str(), colorLegend, colorBackGround, 1, 10);
      
      ////////////////////////////////////////////////////////////////////
      // Axes
      ////////////////////////////////////////////////////////////////////
      
      if(drawAxe)
      {
        const unsigned char colorAxes[] = {0,0,0};
        std::cout << "Draw axes" << std::endl;


        img.draw_point(0, (dimY-1)/3 ,colorAxes);
        img.draw_point(0, 2*(dimY-1)/3 ,colorAxes);
        img.draw_point(dimX-1, (dimY-1)/3 ,colorAxes);
        img.draw_point(dimX-1,  2*(dimY-1)/3 ,colorAxes);
        img.draw_point((dimX-1)/3, 0 ,colorAxes);
        img.draw_point(2*(dimX-1)/3, 0 ,colorAxes);
        img.draw_point((dimX-1)/3, dimY-1 ,colorAxes);
        img.draw_point(2*(dimX-1)/3, dimY-1 ,colorAxes);



/*
        // X
        // 0
        img.draw_point(margeX,  dimY + margeY +0,colorAxes);
        img.draw_point(margeX,  dimY + margeY +1,colorAxes);
        img.draw_point(margeX,  dimY + margeY +2,colorAxes);

        // dimX
        img.draw_point(margeX + dimX-1,  dimY + margeY +0,colorAxes);
        img.draw_point(margeX + dimX-1,  dimY + margeY +1,colorAxes);
        img.draw_point(margeX + dimX-1,  dimY + margeY +2,colorAxes);

        
        // Y
        // 0
        img.draw_point(margeX-1,  dimY + margeY ,colorAxes);
        img.draw_point(margeX-2,  dimY + margeY,colorAxes);
        img.draw_point(margeX-3,  dimY + margeY,colorAxes);
        
        // dimY
//        img.draw_point(margeX + dimX,  dimY + margeY +1,colorAxes);
//        img.draw_point(margeX + dimX,  dimY + margeY +2,colorAxes);
*/      

      }
      
      if(display)
      {
        img.display(title.c_str());
      }
      img.save_jpeg(outImgName.c_str());
      file.close();
      
      std::cout << std::endl << "done." << std::endl;
      std::cout << "Files generated : " << outImgName << std::endl;
      std::cout << "                  " << fileLogName << std::endl;
    }
    else
    {
      std::cerr << "Error while opening the file." << std::endl;
    }
  }
  else
  {
    std::cerr << "Usage : hyplot ./filename" << std::endl;
  }
  return 0;
}
