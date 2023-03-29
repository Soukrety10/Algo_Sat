#include <iostream>
#include "stdafx.h"
#include <iomanip>
#include <stdio.h>
#include <chrono>
#include <ctime>
#include "coreLib.h"
#include "cJulian.h"
#include <vector>
#include "orbitLib.h"

using namespace Zeptomoby::OrbitTools;

using namespace std;


void PrintPosVel(const cSatellite& sat);

bool isVisible(double a, double e, double R, double theta, double min_elevation, double i) {
   if (theta < i){
      return true;
   } 
   else{
   theta = deg2rad(theta); 
   i = deg2rad(i); 
   min_elevation = deg2rad(min_elevation);

   double inter = (R*cos(theta) - a*cos(i)) / tan(theta + min_elevation);
   double thresh = inter + R*sin(theta);
   return a*sin(i) >= thresh;}
 
     
}

//Give The Category of the code
std::string determineVisibilityCategory(double max_elevation) {
    if (max_elevation < 10.0)
        return "Not visible";
    else if (max_elevation < 30.0)
        return "Marginal";
    else if (max_elevation < 60.0)
        return "Good";
    else
        return "Excellent";
}


// Direction from azimuth to NSEW notation
std::string azimuthToNSEW(double azimuth) {
    std::string direction;

    if (azimuth >= 337.5 || azimuth < 22.5) {
        direction = "N";
    } else if (azimuth >= 22.5 && azimuth < 67.5) {
        direction = "NE";
    } else if (azimuth >= 67.5 && azimuth < 112.5) {
        direction = "E";
    } else if (azimuth >= 112.5 && azimuth < 157.5) {
        direction = "SE";
    } else if (azimuth >= 157.5 && azimuth < 202.5) {
        direction = "S";
    } else if (azimuth >= 202.5 && azimuth < 247.5) {
        direction = "SW";
    } else if (azimuth >= 247.5 && azimuth < 292.5) {
        direction = "W";
    } else if (azimuth >= 292.5 && azimuth < 337.5) {
        direction = "NW";
    }

    return direction;
}



// Find when the satelitte will be above the minimum elevation 
std::vector<std::tuple<cJulian, cJulian, std::string, double, double, double, double, double, double, std::string, std::string, std::string>> findNextPasses(const cSatellite& sat, const cSite& site, const cJulian& start_jd, int max_minutes, double minElevation) {
    std::vector<std::tuple<cJulian, cJulian, std::string, double, double, double, double, double, double, std::string, std::string, std::string>> passes;
   

   double max_elevation = 0.0; 
   double max_azimuth = 0.0; 
   bool pass_found = false;
   double elevation_start = 0.0;
   double elevation_end = 0.0; 
   double azimuth_start = 0.0; 
   double azimuth_end = 0.0;
   std::string azimuth_start_direction; 
   std::string azimuth_end_direction;
   std::string azimuth_max_direction;
   cJulian pass_start;

    for (int i = 0; i < max_minutes; ++i) {
        cJulian curr_jd(start_jd);
        curr_jd.AddMin(i);
        cEciTime eci = sat.PositionEci(curr_jd);
        cTopo topoLook = site.GetLookAngle(eci);

        double elevation = topoLook.ElevationDeg();
        double azimuth = topoLook.AzimuthDeg();

        if (elevation >= minElevation ) {

            if (!pass_found) {
                elevation_start  = elevation; 
                azimuth_start = azimuth;
                azimuth_start_direction = azimuthToNSEW(azimuth);
                pass_start = curr_jd;
                pass_found = true;
            }
             if(elevation >= max_elevation){
                   max_elevation = elevation;
                   max_azimuth = azimuth; 
                   azimuth_max_direction = azimuthToNSEW(azimuth);

               }
        } else if (pass_found) {
            cJulian pass_end = curr_jd;
            std::string visibility_category = determineVisibilityCategory(max_elevation);
            passes.push_back(std::make_tuple(pass_start, pass_end, visibility_category, max_elevation, max_azimuth, elevation_start, azimuth_start, elevation_end, azimuth_end, azimuth_start_direction, azimuth_max_direction, azimuth_end_direction));
            pass_found = false;
            max_elevation = 0.0; 
            max_azimuth = 0.0;
            elevation_end  = elevation; 
            azimuth_end = azimuth;
            azimuth_end_direction = azimuthToNSEW(azimuth);  
            
        }
    }

    return passes;
}



// Main function
int main(int argc, char **argv) {
    // Test SGP4 TLE data
    string str1 = "SGP4 Test";
    string str2 = "1 25544U 98067A   23088.63976610  .00039281  00000-0  69794-3 0  9990";
    string str3 = "2 25544  51.6427   2.6677 0005712 146.1232 286.9579 15.49778825389450";

      // Extract mean motion from line 2
    double mean_motion = stod(str2.substr(52, 11));

    // Extract eccentricity from line 2
    double eccentricity = stod("0." + str2.substr(26, 7));

    // Extract inclination from line 3
    double inclination = stod(str3.substr(8, 8));

    

    // Create a TLE object using the data above
    cTle tleSGP4(str1, str2, str3);

    // Create a satellite object from the TLE object
    cSatellite satSGP4(tleSGP4);

    // Print the position and velocity information of the satellite
    PrintPosVel(satSGP4);

    printf("Example output:\n");

   // 90 minutes more than the epock time (it's just a random position)
    cEciTime eciSGP4 = satSGP4.PositionEci(90.0);

    // Coordinates of the observator
    double lon = 48.7107;  
    double lat = 2.1673;
    cSite siteEquator(lon, lat, 0); // N,  W, 0 km altitude

   

    cTopo topoLook = siteEquator.GetLookAngle(eciSGP4);

    // Print out the results.
    printf("AZ: %.3f  EL: %.3f\n", topoLook.AzimuthDeg(), topoLook.ElevationDeg());

    std::time_t t = std::time(0);
    
    std::tm *now = std::gmtime(&t);

 
    std::time_t now_time_t = mktime(now); // Convert std::tm to time_t

    // Get the Julian date
    Zeptomoby::OrbitTools::cJulian start_jd(now_time_t);

    
    int max_minutes = 2*24 * 60; // Find the next satellite passes within the next 2 days

   const double G = 6.6743e-11; // Gravitational constant
   const double M = 5.9722e24; // Mass of the Earth in kg
   const double R = 6371.0e3; // Mean radius of the Earth in meters


   double n = mean_motion * 2 * PI /86400.0;
   double a = pow(G*M/pow(n, 2), 1/3); 

    double min_elevation = 10; 

   if (! isVisible( a,  eccentricity,  R, lat, min_elevation, inclination)){
    
      std::cout << " Satelite is not visible from this location"  <<  std::endl;
      
      }

   else{ auto passes = findNextPasses(satSGP4, siteEquator, start_jd, max_minutes, min_elevation); // The last parameter is the minimum elevation degree

    // Display the passes
    for (const auto &[next_pass_start, next_pass_end, visibility_category, max_elevation, azimuth, elevation_start, azimuth_start, elevation_end, azimuth_end, azimuth_start_direction, azimuth_max_direction, azimuth_end_direction ] : passes) {
        time_t start_time = next_pass_start.ToTime() + (2 * 60 * 60); // To make it UTC + 2
        time_t end_time = next_pass_end.ToTime() + (2 * 60 * 60);

      
        std::tm start_tm;
        gmtime_s(&start_tm, &start_time);

        std::tm end_tm;
        gmtime_s(&end_tm, &end_time);

        std::cout << "Next pass start: " << std::put_time(&start_tm, "%Y-%m-%d %H:%M:%S") << std::endl;
        std::cout << "Next pass start elevation: " << elevation_start << std::endl;
        std::cout << "Next pass start azimuth: " << azimuth_start << std::endl;
        std::cout << "Next pass start azimuth diretion: " << azimuth_start_direction << std::endl;
        std::cout << "Next pass end: " << std::put_time(&end_tm, "%Y-%m-%d %H:%M:%S") << std::endl;
        std::cout << "Next pass end elevation: " << elevation_end << std::endl;
        std::cout << "Next pass end azimuth: " << azimuth_end << std::endl;
         std::cout << "Next pass end azimuth diretion: " << azimuth_end_direction << std::endl;
        std::cout << "Visibility category: " << visibility_category << std::endl;
        std::cout << "Maximum elevation: " << max_elevation << std::endl;
        std::cout << "Maximum Azimuth: " << azimuth << std::endl;
         std::cout << "Next pass maximum azimuth diretion: " << azimuth_max_direction << std::endl;
        std::cout << "------------------------------------------" << std::endl;
    }}

    return 0;
}


// Helper function to output position and velocity information
void PrintPosVel(const cSatellite& sat)
{
   vector<cEci> vecPos;

   // Calculate the position and velocity of the satellite for various times.
   // mpe = "minutes past epoch"
   for (int mpe = 0; mpe <= (360 * 4); mpe += 360)
   {
      // Get the position of the satellite at time "mpe"
      cEciTime eci = sat.PositionEci(mpe);
    
      // Push the coordinates object onto the end of the vector.
      vecPos.push_back(eci);
   }

   // Print TLE data
   printf("%s\n",   sat.Name().c_str());
   printf("%s\n",   sat.Orbit().TleLine1().c_str());
   printf("%s\n\n", sat.Orbit().TleLine2().c_str());

   // Header
   printf("  TSINCE            X                Y                Z\n\n");

   // Iterate over each of the ECI position objects pushed onto the
   // position vector, above, printing the ECI position information
   // as we go.
   for (unsigned int i = 0; i < vecPos.size(); i++)
   {
      printf("%8d.00  %16.8f %16.8f %16.8f\n",
               i * 360,
               vecPos[i].Position().m_x,
               vecPos[i].Position().m_y,
               vecPos[i].Position().m_z);
   }

   printf("\n                    XDOT             YDOT             ZDOT\n\n");

   // Iterate over each of the ECI position objects in the position
   // vector again, but this time print the velocity information.
   for (unsigned int i = 0; i < vecPos.size(); i++)
   {
      printf("             %16.8f %16.8f %16.8f\n",
             vecPos[i].Velocity().m_x,
             vecPos[i].Velocity().m_y,
             vecPos[i].Velocity().m_z);
   }

   printf("\n");

  
}

