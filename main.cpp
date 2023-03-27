
#include "stdafx.h"

#include <stdio.h>
#include <chrono>
#include <ctime>

#include "coreLib.h"


#include "orbitLib.h"

// Forward declaration of helper function; see below
void PrintPosVel(const cSatellite& sat);
// Helper function to find the pass time when the satellite is at least minElevation degrees above the horizon and within a specified azimuth range
double FindPassTime(const cSatellite& sat, const cSite& site, double minElevation, double minAzimuth, double maxAzimuth)
{
   double timeStep = 1.0; // in minutes

   // Get the current time in seconds since the UNIX epoch (Jan 1, 1970)
   auto now = std::chrono::system_clock::now();
   std::time_t currentTimeSec = std::chrono::system_clock::to_time_t(now);
   double currentTime = std::difftime(currentTimeSec, 0) / 60.0;

    return currentTimeSec;


   // Iterate over time steps to find when the satellite is above minElevation and within the specified azimuth range
   /*while (currentTime <= (360 * 4)) // search for 4 orbits (roughly 24 hours)
   {
      // Get the position of the satellite at the current time
      cEciTime eci = sat.PositionEci(currentTime);

      // Get the look angle from the site to the satellite
      cTopo topoLook = site.GetLookAngle(eci);

      // Check if the elevation is greater than or equal to the minimum elevation and if the azimuth is within the specified range
      if (topoLook.ElevationDeg() >= minElevation && topoLook.AzimuthDeg() >= minAzimuth && topoLook.AzimuthDeg() <= maxAzimuth)
      {
          return currentTime;
      }

      // Increment the current time by the time step
      currentTime += timeStep;
   }

   // Return -1 if no suitable pass time is found
   return -1.0;*/
}


//////////////////////////////////////////////////////////////////////////////
int main(int  argc , char **  argv)
{
   // Test SGP4 TLE data
   string str1 = "SGP4 Test";
   string str2 = "1 25544U 98067A   23082.39717484  .00016571  00000+0  30373-3 0  9996";
   string str3 = "2 25544  51.6420  33.5716 0006039 122.8805  22.8850 15.49417071388484";

   // Create a TLE object using the data above
   cTle tleSGP4(str1, str2, str3);

   // Create a satellite object from the TLE object
   cSatellite satSGP4(tleSGP4);

   // Print the position and velocity information of the satellite
   PrintPosVel(satSGP4);

   

  

   printf("Example output:\n");

   cEciTime eciSGP4 = satSGP4.PositionEci(90.0);

   
   cSite siteEquator(-26.5816902, 18.1365094, 0); // 0.00 N, 100.00 W, 0 km altitude

  
   cTopo topoLook = siteEquator.GetLookAngle(eciSGP4);

   // Print out the results.
   printf("AZ: %.3f  EL: %.3f\n", 
          topoLook.AzimuthDeg(),
          topoLook.ElevationDeg());

   double minAzimuth = 0.0;   // example minimum azimuth angle
   double maxAzimuth = 180.0; // example maximum azimuth angle

    // Find the pass time when the satellite is at least 15 degrees above the horizon and within the specified azimuth range
    double passTime = FindPassTime(satSGP4, siteEquator, 15.0, minAzimuth, maxAzimuth);

    printf("Pass time with elevation >= 15.0 degrees and azimuth between %.2f and %.2f degrees: %.2f minutes past epoch\n", minAzimuth, maxAzimuth, passTime);
}

/////////////////////////////////////////////////////////////////////////////
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

