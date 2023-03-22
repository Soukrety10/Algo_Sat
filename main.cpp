#include <iostream>
#include <cmath>


using namespace std;

const double PI = 3.14159265358979323846;
const double GM = 398600.4418; // Constante de gravitation de la terre
const double min_elevation= 15.0;  // Minimum elevation angle in degrees
const double RT = 6371.0; // Rayon de la Terre

// degrees => radians
double deg2rad(double deg) {
return deg * PI / 180.0;
}

double rad2deg(double rad) {
    return rad * 180.0 / PI;
}

double mod(double a, double b) {
    return a - b * floor(a / b);
}

void getPosition(double inclination, double eccentricity, double meanMotion, double raaa, double argumentPerigee, double meanAnomaly, double epochTime, double currentTime, double& x, double& y, double& z) {
    
    // Convertion
    inclination = deg2rad(inclination);
    raaa = deg2rad(raaa);
    argumentPerigee  = deg2rad(argumentPerigee);
    meanAnomaly = deg2rad(meanAnomaly);
    double n = meanMotion * 2 * PI / 86400.0; // Vitesse angulaire en rad/s

    double a = pow((GM / (n*n)), 1.0 / 3.0);   // 3éme loi de Kepler
    

    // Anomalie Actuel (angle du satelite dans sa position actuel dans son orbite avec son éventuel position dans un orbite circulaire parfait )

    double epochAnomaly = meanAnomaly;
    double deltaT = currentTime - epochTime;

    double currentAnomaly = epochAnomaly + (n * deltaT); // n*deltaT est l'angle balayé par le satelitte dans son mouvement

    // Distance entre le Satelitte et le centre de la terre
    double r = a * (1 - eccentricity * cos(currentAnomaly));
    double trueAnomaly = currentAnomaly + argumentPerigee ;

    // Calcul des Coordonnées dans le référentiel terestre (projection de lambert)
    x = r * (cos(raaa) * cos(trueAnomaly) - sin(raaa) * sin(trueAnomaly) * cos(inclination));
    y = r * (sin(raaa) * cos(trueAnomaly) + cos(raaa) * sin(trueAnomaly) * cos(inclination));
    z = r * (sin(trueAnomaly) * sin(inclination));
}

double calculateElevation(double obs_lat, double obs_lon, double obs_alt, double x, double y, double z) {
    double obs_x = (RT + obs_alt) * cos(deg2rad(obs_lat)) * cos(deg2rad(obs_lon));
    double obs_y = (RT + obs_alt) * cos(deg2rad(obs_lat)) * sin(deg2rad(obs_lon));
    double obs_z = (RT + obs_alt) * sin(deg2rad(obs_lat));

    double range_x = x - obs_x;
    double range_y = y - obs_y;
    double range_z = z - obs_z;

    double range = sqrt(range_x * range_x + range_y * range_y + range_z * range_z);

    double elevation_rad = asin(range_z / range);
    return rad2deg(elevation_rad);
}

time_t findNextPass(double inclination, double eccentricity, double meanMotion, double raaa, double argumentPerigee, double meanAnomaly, double epochTime, double obs_lat, double obs_lon, double obs_alt, double min_elevation) {
    double currentTime = epochTime;
    double step = 600.0;  // Time step in seconds

    while (true) {
        double x, y, z;
        getPosition(inclination, eccentricity, meanMotion, raaa, argumentPerigee, meanAnomaly, epochTime, currentTime, x, y, z);
        double elevation = calculateElevation(obs_lat, obs_lon, obs_alt, x, y, z);

        if (elevation >= min_elevation) {
            break;
        }

        currentTime += step;
    }

    return static_cast<time_t>(currentTime);
}

int main()
{

    // Boucle sur un fichier 
    //TLE data
    string line1 = "1 25544U 98067A   19079.78876710  .00000786  00000-0  20960-4 0  9992";
    string line2 = "2 25544  51.6426 357.7587 0004698  15.3385  88.5065 15.53281465158156";

    // définition des différents paramétres
    string name;
    double inclination, raan, eccentricity, argumentPerigee, meanAnomaly, meanMotion, epochTime  ;

    
    name = line1.substr(2, 7); // nom du satelite 
    inclination = stod(line2.substr(8, 8)); // Angle Inclinaison  
    raan = stod(line2.substr(17, 8));  // Right ascension
    eccentricity = stod("." + line2.substr(26, 7)); // Eccentricity
    argumentPerigee = stod(line2.substr(34, 8)); // Angle Perigee 
    meanAnomaly = stod(line2.substr(43, 8));// Angle Anomalie Moyenne 
    meanMotion = stod(line2.substr(52, 11));  // Vitesse angulaire 
    epochTime = stod(line1.substr(17, 15)); // Temps de la dérnier validation du position de satelite

    

    

    double min_elevation = 15.0;
     

    // Location des observers ¨Paris
    double obs_lat = 48.8588897;  
    double obs_lon = 2.320041; 
    double obs_alt = 0.0;
    
    double x, y, z; 
    time_t nextPassTime = findNextPass(inclination, eccentricity, meanMotion, raan, argumentPerigee, meanAnomaly, epochTime, obs_lat, obs_lon, obs_alt, min_elevation);

    // Convert the next pass time to a human-readable format
    char timeBuffer[26];
    ctime_s(timeBuffer, sizeof timeBuffer, &nextPassTime);
    string nextPassTimeString = string(timeBuffer);

    // Print the time of the next pass
    cout << "The next pass of satellite " << name << " over the observer's location will be at: " << nextPassTimeString << endl;

    return 0;
}
