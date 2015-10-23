//
//  main.cpp
//  CG_hw5
//
//  Created by Shangqi Wu on 14/11/13.
//  Copyright (c) 2014 Shangqi Wu. All rights reserved.
//

#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
#include <limits>
#include <algorithm>

using namespace std;

// Claiming classes and structs will be used
class Point;
class Vertex;
class Line;
class Polygon;
struct Face{
    int p1;
    int p2;
    int p3;
};
struct scanline{
    int x;
    double buffered_z;
};
struct pixel{
    char color;
    double buffered_z;
};

// Claiming global variants will be used
int j=0, k=0, o=500, p=500;
double x=0, y=0, z=1, X=0, Y=0, Z=0, q=0, r=0, w=-1, Q=0, R=1, W=0, u=-0.8, v=-0.8, U=0.8, V=0.8, F=0.6, B=-0.6;
string f[3] = {"", "", ""}, output = "out.xpm";
bool P = false; // If it is true, this program uses parallel projection, otherwise it uses perspective projection.
//bool P = true; // For debugging use, testing parallel projection feature.
double front_clip;
vector<vector<pixel> > xpmpixel;

// Claiming functions will be used
string setheader();
string setend();
void help();
int rnd(double arg);
void outfile(string output);
void optana(int argc, char * const argv[]);
void readfile(string input, vector<Vertex> *vecver, vector<Face> *vecf);
vector<Vertex> projection(vector<Vertex> vecver);
bool scancmp(scanline scan_in1, scanline scan_in2);
Polygon shclip(Polygon argp, int a, int b, int c, int d);

//--------------------------------------------------------------------------------------------------
class Point {
private: // Stores x,y coordinates of x, and the computed z
    int px;
    int py;
    double z_buffer;
public:
    void set(double argx, double argy) {
        px = argx;
        py = argy;
    }
    void zbuffer(double z) {
        z_buffer = z;
    }
    int getx() {
        return px;
    }
    int gety() {
        return py;
    }
    double getz() {
        return z_buffer;
    }
    void projtoviewport(double ver_x, double ver_y){ // Translate vertex coordinates to viewport.
        if (P == false) {
            double zproj = abs(z/(B-z));
            double sx = (double)(o-j) / (double)(2*zproj);
            double sy = (double)(p-k) / (double)(2*zproj);
            double matrix[3][3] = {1, 0, zproj, 0, 1, zproj, 0, 0, 1};
            double pt[3] = {0, 0, 0};
            for (int i = 0; i < 3; i++) { // Translate to origin of the world coordinates first.
                pt[i] += matrix[i][0] * ver_x;
                pt[i] += matrix[i][1] * ver_y;
                pt[i] += matrix[i][2];
            }
            double matrix_1[3][3] = {sx, 0, 0, 0, sy, 0, 0, 0, 1};
            double pt_1[3] = {0, 0, 0};
            for (int i = 0; i < 3  ; i++) { // Scale coordinates to fit viewport.
                pt_1[i] += matrix_1[i][0] * pt[0];
                pt_1[i] += matrix_1[i][1] * pt[1];
                pt_1[i] += matrix_1[i][2];
            }
            double matrix_2[3][3] = {1, 0, (double)j, 0, 1, (double)k, 0, 0, 1};
            double pt_2[3] = {0, 0, 0};
            for (int i = 0; i < 3; i++) { // Translate point to viewport.
                pt_2[i] += matrix_2[i][0] * pt_1[0];
                pt_2[i] += matrix_2[i][1] * pt_1[1];
                pt_2[i] += matrix_2[i][2];
            }
            px = rnd(pt_2[0]); py = rnd(pt_2[1]);
        }else{ // This part is for parallel projection, first part translates to origin by 2, the rest are the same as perspective projection.
            double sx = (double)(o-j) / (double)(2);
            double sy = (double)(p-k) / (double)(2);
            double matrix[3][3] = {1, 0, 1, 0, 1, 1, 0, 0, 1};
            double pt[3] = {0, 0, 0};
            for (int i = 0; i < 3; i++) {
                pt[i] += matrix[i][0] * ver_x;
                pt[i] += matrix[i][1] * ver_y;
                pt[i] += matrix[i][2];
            }
            double matrix_1[3][3] = {sx, 0, 0, 0, sy, 0, 0, 0, 1};
            double pt_1[3] = {0, 0, 0};
            for (int i = 0; i < 3  ; i++) {
                pt_1[i] += matrix_1[i][0] * pt[0];
                pt_1[i] += matrix_1[i][1] * pt[1];
                pt_1[i] += matrix_1[i][2];
            }
            double matrix_2[3][3] = {1, 0, (double)j, 0, 1, (double)k, 0, 0, 1};
            double pt_2[3] = {0, 0, 0};
            for (int i = 0; i < 3; i++) {
                pt_2[i] += matrix_2[i][0] * pt_1[0];
                pt_2[i] += matrix_2[i][1] * pt_1[1];
                pt_2[i] += matrix_2[i][2];
            }
            px = rnd(pt_2[0]); py = rnd(pt_2[1]);
        }
    }
    bool operator==(Point a) {
        return (px==a.getx() && py==a.gety());
    }
    bool operator!=(Point a) {
        return (px!=a.getx() || py!=a.gety());
    }
};

//-------------------------------------------------------------------------------------------------
class Vertex{
private:
    double verx;
    double very;
    double verz;
    int status; // Indicating positon of the vertex by Cohen-Sutherland code.
public:
    void set(double a, double b, double c){
        verx = a;
        very = b;
        verz = c;
    }
    int get_status(){
        return status;
    }
    void trivial_test(){
        status = 0b00000000;
        if (P == true) { // Trivial test for parallel projection.
            if (very > 1) { // Using 6-bit code for Cohen-Sutherland trivial test
                status |= 0b000001;
            }else if (very < -1){
                status |= 0b000010;
            }
            if (verx > 1) {
                status |= 0b000100;
            }else if (verx < -1){
                status |= 0b001000;
            }
            if (verz < -1) {
                status |= 0b010000;
            }else if (verz > 0){
                status |= 0b100000;
            }
        }else{
            if (very > -verz) {
                status |= 0b000001;
            }else if (very < verz){
                status |= 0b000010;
            }
            if (verx > -verz) {
                status |= 0b000100;
            }else if (verx < verz){
                status |= 0b001000;
            }
            if (verz < -1) {
                status |= 0b010000;
            }else if (verz > front_clip){
                status |= 0b100000;
            }
        }
    }
    void matrix_multiply(double a[4][4]){ // To implement the multiplication with normalizing matrix.
        double pt[4] = {0, 0, 0, 0};
        double temp = 0;
        for (int i = 0; i < 4; i++) {
            temp = temp + (a[i][0] * verx + a[i][1] * very);
            temp = temp + (a[i][2] * verz + a[i][3]);
            pt[i] = temp;
            temp = 0;
        }
        verx = pt[0]; very = pt[1]; verz = pt[2];
    }
    Point toProjection(){
        Point a;
        if (P == true) { // If implementing parallel projection, directly delete z coordinates.
            a.projtoviewport(verx, very);
            a.zbuffer(verz);
        }else{ // Implementing perspective projection to the desired plane.
            double zproj = z / (B-z);
            if (verz <= front_clip) {
                double xtmp = verx * zproj / verz;
                double ytmp = very * zproj / verz;
                a.projtoviewport(xtmp, ytmp);
                a.zbuffer(verz);
            }else if (front_clip < verz) { // Still generate a boundry point within the world window, to ensure continuity of the image.
                double xtmp = verx * zproj / front_clip;
                double ytmp = very * zproj / front_clip;
                a.projtoviewport(xtmp, ytmp);
                a.zbuffer(verz);
            }
        }
        return a;
    }
};

//--------------------------------------------------------------------------------------------------
class Line { // present the line by: y = slope*x + d
private:
    Point start;
    Point end;
    int xmax;
    int xmaxy;
    int xmin;
    int xminy;
    int ymax;
    int ymin;
    bool sd_exist;
    double slope;
    double d1;
    bool in;
    double xminz;
    double xmaxz;
    void cal() {
        if (start.getx() != end.getx()) {
            if (start.getx() > end.getx()) {
                xmax = start.getx();
                xmaxy = start.gety();
                xmin = end.getx();
                xminy = end.gety();
                xminz = end.getz();
                xmaxz = start.getz();
            }else {
                xmax = end.getx();
                xmaxy = end.gety();
                xmin = start.getx();
                xminy = start.gety();
                xminz = start.getz();
                xmaxz = end.getz();
            }
            if (start.gety()>=end.gety()) { // set ymax and ymin
                ymax = start.gety();
                ymin = end.gety();
            }else {
                ymax = end.gety();
                ymin = start.gety();
            }
            slope = (double)(start.gety() - end.gety()) / (double)(start.getx() - end.getx());
            d1 = start.gety() - (slope * start.getx());
            sd_exist = true;
        }else {
            if (start.gety()>=end.gety()) {
                ymax = start.gety();
                ymin = end.gety();
                xmaxy = start.gety();
                xminy = end.gety();
                xminz = end.getz();
                xmaxz = start.getz();
            }else {
                ymax = end.gety();
                ymin = start.gety();
                xmaxy = end.gety();
                xminy = start.gety();
                xminz = start.getz();
                xmaxz = end.getz();
            }
            xmax = start.getx();
            xmin = end.getx();
            slope = numeric_limits<double>::max();
            d1 = numeric_limits<double>::max();
            sd_exist = false;
        }
    }
public:
    void setp(Point a, Point b) {
        start = a;
        end = b;
        cal();
    }
    bool sd() {
        return sd_exist;
    }
    Point gets() {
        return start;
    }
    Point gete() {
        return end;
    }
    int cal_inter_wlr(int x) { // return the intersaction point of x=x
        if (sd_exist) {
            double y = slope*(double)x + d1;
            return rnd(y);
        }else {
            return numeric_limits<int>::max();
        }
    }
    int cal_inter_wtb(int y) { // return the intersaction point of y=y
        if (sd_exist && slope!=0) {
            double x = ((double)y - d1) / slope;
            return rnd(x);
        }else if (slope == 0) {
            return numeric_limits<int>::max();
        }else {
            return start.getx();
        }
    }
    void setin(bool ifin) {
        in = ifin;
    }
    double getslope() {
        return slope;
    }
    bool getin() {
        return in;
    }
    int getxmin() {
        return xmin;
    }
    int getxmax() {
        return xmax;
    }
    int getymin() {
        return ymin;
    }
    int getymax() {
        return ymax;
    }
    int getxmaxy() {
        return xmaxy;
    }
    int getxminy() {
        return xminy;
    }
    double getzleft() {
        return xminz;
    }
    double getzright() {
        return xmaxz;
    }
    vector<vector<scanline> > scanline_register(vector<vector<scanline> > scan_input) { // Calculate all the points of the line.
        scanline new_scan;
        if (sd_exist == true) {
            // Fllowing codes are of Bresenham Algorithm, presented in L-02_Lines.pdf. Codes are modifiied for this cpp file.
            int dx, dy, D, x_br, y_br;
            dx = xmax - xmin;
            dy = ymax - ymin;
            if (0<slope && slope<1) {
                D = 2*dy - dx;
                y_br = ymin;
                for (x_br = xmin; x_br <= xmax; x_br++) {
                    if (D <= 0) {
                        D += 2*dy;
                    }else {
                        D += 2*(dy - dx);
                        if (y_br != ymax) { // Use information of edge, add x and z coordinates on each y as filling boundary.
                            new_scan.x = x_br;
                            new_scan.buffered_z = xminz - (xminz-xmaxz)*(xmin-x_br)/(xmin-xmax);
                            scan_input[y_br].push_back(new_scan);
                        }
                        y_br++;
                    }
                }
            }else if (slope > 1) {
                D = 2*dx - dy;
                x_br = xmin;
                for (y_br = ymin; y_br <= ymax; y_br++) {
                    if (y_br != ymax) {
                        new_scan.x = x_br;
                        new_scan.buffered_z = xminz - (xminz-xmaxz)*(ymin-y_br)/(ymin-ymax);
                        scan_input[y_br].push_back(new_scan);
                    }
                    if (D <= 0) {
                        D += 2*dx;
                    }else {
                        D += 2*(dx - dy);
                        x_br++;
                    }
                }
            }else if (-1<slope && slope<0) {
                D = 2*dy - dx;
                y_br = ymin;
                for (x_br = xmax; x_br >= xmin; x_br--) {
                    if (D <= 0) {
                        D += 2*dy;
                    }else {
                        if (y_br != ymax) {
                            new_scan.x = x_br;
                            new_scan.buffered_z = xminz - (xminz-xmaxz)*(xmin-x_br)/(xmin-xmax);
                            scan_input[y_br].push_back(new_scan);
                        }
                        D += 2*(dy - dx);
                        y_br++;
                    }
                }
            }else if (slope == 1) {
                y_br = ymin;
                for (x_br = xmin; x_br <= xmax; x_br++) {
                    if (y_br != ymax) {
                        new_scan.x = x_br;
                        new_scan.buffered_z = xminz - (xminz-xmaxz)*(xmin-x_br)/(xmin-xmax);
                        scan_input[y_br].push_back(new_scan);
                    }
                    y_br++;
                }
            }else if (slope == -1) {
                y_br = ymin;
                for (x_br = xmax; x_br >= xmin; x_br--) {
                    y_br++;
                    if (y_br != ymax) {
                        new_scan.x = x_br;
                        new_scan.buffered_z = xminz - (xminz-xmaxz)*(xmin-x_br)/(xmin-xmax);
                        scan_input[y_br].push_back(new_scan);
                    }
                }
            }else if (slope == 0) {
               // By the filling algorithm, just ignore horizontal lines.
            }else { // i.e., slope<-1
                D = 2*dx - abs(dy);
                x_br = xmin;
                for (y_br = ymax; y_br >= ymin; y_br--) {
                    if (y_br != ymax) {
                        new_scan.x = x_br;
                        new_scan.buffered_z = xminz - (xminz-xmaxz)*(ymax-y_br)/(ymax-ymin);
                        scan_input[y_br].push_back(new_scan);
                    }
                    if (D <= 0) {
                        D += 2*dx;
                    }else {
                        D += 2*(dx - dy);
                        x_br++;
                    }
                }
            }
        }else if (sd_exist == false) { // for vertical lines
            int x_br = xmin;
            for (int y_br = ymin; y_br <= ymax; y_br++) {
                if (y_br != ymax) {
                    new_scan.x = x_br;
                    new_scan.buffered_z = xminz - (xminz-xmaxz)*(ymin-y_br)/(ymin-ymax);
                    scan_input[y_br].push_back(new_scan);
                }
            }
        }
        return scan_input;
    }
};

class Polygon{ // A set of all vertices and edges of a polygon, with necessary function.
private:
    vector<Point> vertice;
    vector<Line> edge;
    vector<vector<scanline> > scanreg;
    bool out;
    int xmin;
    int xmax;
    int ymin;
    int ymax;
    void calp() {
        edge.resize(vertice.size()-1);
        xmin = vertice[0].getx();
        xmax = xmin;
        ymin = vertice[0].gety();
        ymax = ymin;
        for (int i = 0; i < vertice.size()-1; i++) {
            edge[i].setp(vertice[i], vertice[i+1]);
            if (vertice[i].getx() < xmin) {
                xmin = vertice[i].getx();
            }
            if (vertice[i].getx() > xmax) {
                xmax = vertice[i].getx();
            }
            if (vertice[i].gety() < ymin) {
                ymin = vertice[i].gety();
            }
            if (vertice[i].gety() > ymax) {
                ymax = vertice[i].gety();
            }
        }
        scanreg.resize(501);
    }
    void call() {
        vertice.resize(edge.size());
        xmax = edge[0].gets().getx();
        xmin = xmax;
        ymax = edge[0].gets().gety();
        ymin = ymax;
        for (int i = 0; i < edge.size(); i++) {
            vertice[i] = edge[i].gets();
            if (vertice[i].getx() < xmin) {
                xmin = vertice[i].getx();
            }
            if (vertice[i].getx() > xmax) {
                xmax = vertice[i].getx();
            }
            if (vertice[i].gety() < ymin) {
                ymin = vertice[i].gety();
            }
            if (vertice[i].gety() > ymax) {
                ymax = vertice[i].gety();
            }
        }
        vertice.push_back(vertice[0]);
        scanreg.resize(501);
    }
public:
    void setp(vector<Point> argp) { // input vertices of a polygon and then calculate its edges
        vertice.clear();
        edge.clear();
        if (argp.size() > 1) {
            if (argp[0] != argp[argp.size()-1]) {
                argp.push_back(argp[0]);
            } // If you choose to set up a polygon by points, the first point and the last must be the same.
            for (int i = 0; i < argp.size(); i++) {
                vertice.push_back(argp[i]);
            }
            calp();
            out = false;
        }else if (argp.size() == 1) {
            vertice.push_back(argp[0]);
            vertice.push_back(argp[0]);
            out = false;
        }else {
            out = true;
        }
    }
    void setl(vector<Line> argl) { // input lines to calculate vertices
        edge.clear();
        vertice.clear();
        if (argl.size() != 0) {
            if (argl[0].gets() != argl[argl.size()-1].gete()) {
                Line buff;
                buff.setp(argl[argl.size()-1].gete(), argl[0].gets());
                argl.push_back(buff);
            } // If you choose to set up a polygon by edges, starting point of the first edge and ending point of the last edge must be the same.
            edge = argl;
            call();
            out = false;
        }else {
            out = true;
        }
    }
    vector<Point> getparr() { // return all vertices
        return vertice;
    }
    vector<Line> getlarr() { // return all lines
        return edge;
    }
    Point getp(int i) { // return the ith point
        return vertice[i];
    }
    Line getl(int i) { // return the ith line
        return edge[i];
    }
    bool getout() {
        return out;
    }
    int getxmin() {
        return xmin;
    }
    int getxmax() {
        return xmax;
    }
    int getymin() {
        return ymin;
    }
    int getymax() {
        return ymax;
    }
    void pre_fill() { // Get scanline information for each y.
        for (int i = 0; i < edge.size(); i++) {
            scanreg = edge[i].scanline_register(scanreg);
        }
    }
    void fillin(int num) { // Fill in polygons by each scanline
        for (int i = ymin; i < ymax; i++) {
            int size = (int)scanreg[i].size();
            sort(scanreg[i].begin(), scanreg[i].end(),scancmp);
            for (int j = 0; j < size; j+=2) {
                for (int k = scanreg[i][j].x+1; k <= scanreg[i][j+1].x; k++) {
                    double z_pixel = scanreg[i][j].buffered_z-(scanreg[i][j].buffered_z-scanreg[i][j+1].buffered_z)*(scanreg[i][j].x-k)/(scanreg[i][j].x-scanreg[i][j+1].x); // Calculate z coordinates for each pixel.
                    if (z_pixel<front_clip && z_pixel>=xpmpixel[i][k].buffered_z) {
                        xpmpixel[i][k].buffered_z = z_pixel; // If new z is within view volume and closer to screen, replace the old z.
                        double z_ratio = (z_pixel+1)/(front_clip+1);
                        //num = 2;
                        int level = (int)(z_ratio*19.6);
                        if (num == 0) {
                            xpmpixel[i][k].color = (char)(60+level);
                        }else if (num == 1) {
                            xpmpixel[i][k].color = (char)(80+level);
                        }else if (num == 2) {
                            xpmpixel[i][k].color = (char)(100+level);
                        }
                    }
                }
            }
        }
    }
};

// Main Function------------------------------------------------------------------------------------
int main(int argc, char * argv[]) {
    // Prepare pixels and input instructions.
    int height = 501;
    int width = 501;
    xpmpixel.resize(height);
    for (int i = 0; i < height; i++) {
        xpmpixel[i].resize(width);
        for (int j = 0; j < width; j++) {
            xpmpixel[i][j].color = '-';
            xpmpixel[i][j].buffered_z = -1;
        }
    }
    // Analyzing input instructions.
    //optana(argc, argv);
    // Prepare front clipping plane first.
    if (P == false) {
        front_clip = (z-F) / (B-z);
    }else {
        front_clip = 0;
    }
    f[0] = "/Users/wushangqi/cghw4/bound-sprellpsd.smf"; // dir only for debug*****************************************
    //f[0] = "/Users/wushangqi/cghw4/bound-bunny_1k.smf";
    //f[1] = "/Users/wushangqi/cghw4/bound-cow.smf";
    //f[2] = "/Users/wushangqi/cghw4/bbound-sprtrd.smf";
    for (int l = 0; l < 3; l++) {
        if (!f[l].empty()) { // Processing image if the f[l] is not empty. f[0] for red, f[1] for green, and f[2] for blue images.
            vector<Vertex> vecver;
            vector<Face> vecf;
            vector<Polygon> vecpl;
            // Read the file.
            readfile(f[l], &vecver, &vecf);
            // Processing 3D Projection.
            vecver = projection(vecver);
            int fsize = (int)vecf.size();
            // Prepare to draw edges of each face.
            for (int i = 0; i < fsize; i++) { // Analysing trivial test result first.
                int s1 = vecver[vecf[i].p1].get_status() & vecver[vecf[i].p2].get_status();
                int s2 = vecver[vecf[i].p2].get_status() & vecver[vecf[i].p3].get_status();
                int s3 = vecver[vecf[i].p3].get_status() & vecver[vecf[i].p1].get_status();
                if (s1==0 || s2==0 || s3==0) { // If the face is within view volume, translate it into 2-D polygon.
                    vector<Point> vecp;
                    vecp.push_back(vecver[vecf[i].p1].toProjection());
                    vecp.push_back(vecver[vecf[i].p2].toProjection());
                    vecp.push_back(vecver[vecf[i].p3].toProjection());
                    Polygon newpoly;
                    newpoly.setp(vecp);
                    vecpl.push_back(newpoly);
                }
            } // Calculating z information and fillin for each polygon.
            int plsize = (int)vecpl.size();
            for (int i = 0; i < plsize; i++) { // Clip polygon with viewport boundary first.
                vecpl[i] = shclip(vecpl[i], j, k, o, p);
                if (vecpl[i].getout() == false) {
                    vecpl[i].pre_fill();
                    vecpl[i].fillin(l);
                }
            }
        }
    }
    //f[1] = "1";
    //f[2] = "2";
    // Write output xpm file.
    output = "/Users/wushangqi/out.xpm"; // dir only for debug******************************************************
    outfile(output);
    // Display output xpm image automatically.
    string shell = "display " + output;
    system(shell.c_str());
    return 0;
}

//---------------------------------------------------------------------------------------------------
vector<Vertex> projection(vector<Vertex> vecver){
    // Prepare to process vertex.
    int versize = (int)vecver.size();
    vector<double> rz(3); // VPN / |VPN|
    double rzsqrt = sqrt(q*q + r*r + w*w);
    rz[0] = q / rzsqrt;
    rz[1] = r / rzsqrt;
    rz[2] = w / rzsqrt;
    vector<double> rx(3); // VUP x Rz first
    rx[0] = R*rz[2] - W*rz[1];
    rx[1] = W*rz[0] - Q*rz[2];
    rx[2] = Q*rz[1] - R*rz[0];
    double rxsqrt = sqrt(rx[0]*rx[0] + rx[1]*rx[1] + rx[2]*rx[2]); // then normalizing it
    rx[0] = rx[0] / rxsqrt;
    rx[1] = rx[1] / rxsqrt;
    rx[2] = rx[2] / rxsqrt;
    vector<double> ry(3); // Implementing Rz x Rx
    ry[0] = rz[1]*rx[2] - rz[2]*rx[1];
    ry[1] = rz[2]*rx[0] - rz[0]*rx[2];
    ry[2] = rz[0]*rx[1] - rz[1]*rx[0];
    double shx = (0.5*(U+u)-x) / z;
    double shy = (0.5*(V+v)-y) / z;
    // Starting to process 3D projection.
    if (P == true) { // Implementing Parallel Projection
        double tpar1 = -(U+u) / 2.0;
        double tpar2 = -(V+v) / 2.0;
        double spar1 = 2.0 / (U-u);
        double spar2 = 2.0 / (V-v);
        double spar3 = 1.0 / (F-B);
        double R[4][4] = {{rx[0], rx[1], rx[2], 0}, {ry[0], ry[1], ry[2], 0}, {rz[0], rz[1], rz[2], 0}, {0, 0, 0, 1}};
        double Tvrp[4][4] = {{1, 0, 0, -X}, {0, 1, 0, -Y}, {0, 0, 1, -Z}, {0, 0, 0, 1}};
        double SHpar[4][4] = {{1, 0, shx, 0}, {0, 1, shy, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
        double Spar[4][4] = {{spar1, 0, 0, 0}, {0, spar2, 0, 0}, {0, 0, spar3, 0}, {0, 0, 0, 1}};
        double Tpar[4][4] = {{1, 0, 0, tpar1}, {0, 1, 0, tpar2}, {0, 0, 1, -F}, {0, 0, 0, 1}};
        double tmp1[4][4] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        double tmp2[4][4] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        double tmp3[4][4] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        double Npar[4][4] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        double temp = 0;
        for (int i = 0; i < 4; i++) { // Pipeline Step 1: Translate VRP to the origin.
            for (int j = 0; j < 4; j++) { // Pipeline Step 2: Rotate so VPN becomes z, VUP becomes y, and u becomes x.
                temp = temp + (R[i][0] * Tvrp[0][j] + R[i][1] * Tvrp[1][j]);
                temp = temp + (R[i][2] * Tvrp[2][j] + R[i][3] * Tvrp[3][j]);
                tmp1[i][j] = temp;
                temp = 0;
            }
        }
        for (int i = 0; i < 4; i++) { // Pipeline Step 3: Shear to make direction of the projection become parallel to z.
            for (int j = 0; j < 4; j++) {
                temp = temp + (SHpar[i][0] * tmp1[0][j] + SHpar[i][1] * tmp1[1][j]);
                temp = temp + (SHpar[i][2] * tmp1[2][j] + SHpar[i][3] * tmp1[3][j]);
                tmp2[i][j] = temp;
                temp = 0;
            }
        }
        for (int i = 0; i < 4; i++) { // Pipeline Step 4: Translate and scale into a canonical view volume.
            for (int j = 0; j < 4; j++) {
                temp = temp + (Tpar[i][0] * tmp2[0][j] + Tpar[i][1] * tmp2[1][j]);
                temp = temp + (Tpar[i][2] * tmp2[2][j] + Tpar[i][3] * tmp2[3][j]);
                tmp3[i][j] = temp;
                temp = 0;
            }
        }
        for (int i = 0; i < 4; i++) {  // Implementing last scaling step to generate normalizing matrix.
            for (int j = 0; j < 4; j++) {
                temp = temp + (Spar[i][0] * tmp3[0][j] + Spar[i][1] * tmp3[1][j]);
                temp = temp + (Spar[i][2] * tmp3[2][j] + Spar[i][3] * tmp3[3][j]);
                Npar[i][j] = temp;
                temp = 0;
            }
        }
        for (int i = 0; i < versize; i++) { // Using this matrix to process coordinates of all vertex.
            vecver[i].matrix_multiply(Npar);
            vecver[i].trivial_test(); // Implemeting simple view volume rejection test for each vertex.
        }
    }else{ // Implementing Perspective Projection
        double sperx = 2.0*(-z)/((U-u)*(-z+B));
        double spery = 2.0*(-z)/((V-v)*(-z+B));
        double sperz = -1.0/(-z+B);
        double R[4][4] = {{rx[0], rx[1], rx[2], 0}, {ry[0], ry[1], ry[2], 0}, {rz[0], rz[1], rz[2], 0}, {0, 0, 0, 1}};
        double Tvrp[4][4] = {{1, 0, 0, -X}, {0, 1, 0, -Y}, {0, 0, 1, -Z}, {0, 0, 0, 1}};
        double SHper[4][4] = {{1, 0, shx, 0}, {0, 1, shy, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
        double Sper[4][4] = {{sperx, 0, 0, 0}, {0, spery, 0, 0}, {0, 0, sperz, 0}, {0, 0, 0, 1}};
        double Tprp[4][4] = {{1, 0, 0, -x}, {0, 1, 0, -y}, {0, 0, 1, -z}, {0, 0, 0, 1}};
        double tmp1[4][4] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        double tmp2[4][4] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        double tmp3[4][4] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        double Nper[4][4] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        double temp = 0;
        for (int i = 0; i < 4; i++) { // Pipeline Step 1: Translate VRP to the origin.
            for (int j = 0; j < 4; j++) {
                temp = temp + (R[i][0] * Tvrp[0][j] + R[i][1] * Tvrp[1][j]);
                temp = temp + (R[i][2] * Tvrp[2][j] + R[i][3] * Tvrp[3][j]);
                tmp1[i][j] = temp;
                temp = 0;
            }
        }
        for (int i = 0; i < 4; i++) { // Pipeline Step 2: Rotate so VPN becomes z, VUP becomes y, and u becomes x.
            for (int j = 0; j < 4; j++) {
                temp = temp + (Tprp[i][0] * tmp1[0][j] + Tprp[i][1] * tmp1[1][j]);
                temp = temp + (Tprp[i][2] * tmp1[2][j] + Tprp[i][3] * tmp1[3][j]);
                tmp2[i][j] = temp;
                temp = 0;
            }
        }
        for (int i = 0; i < 4; i++) { // Pipeline Step 3: Translate COP to origin.
            for (int j = 0; j < 4; j++) {
                temp = temp + (SHper[i][0] * tmp2[0][j] + SHper[i][1] * tmp2[1][j]);
                temp = temp + (SHper[i][2] * tmp2[2][j] + SHper[i][3] * tmp2[3][j]);
                tmp3[i][j] = temp;
                temp = 0;
            }
        }
        for (int i = 0; i < 4; i++) { // Pipeline Step 5: Scale into a canonical view volume for clippiing.
            for (int j = 0; j < 4; j++) {
                temp = temp + (Sper[i][0] * tmp3[0][j] + Sper[i][1] * tmp3[1][j]);
                temp = temp + (Sper[i][2] * tmp3[2][j] + Sper[i][3] * tmp3[3][j]);
                Nper[i][j] = temp;
                temp = 0;
            }
        } // Generated the normalizing matrix for all coordinates calculation of every vertex.
        for (int i = 0; i < versize; i++) { // Using this matrix to process coordinates of all vertex.
            vecver[i].matrix_multiply(Nper);
            vecver[i].trivial_test(); // Implementing simple view volume rejection test.
        }
    }
    return vecver;
}

//-------------------------------------------------------------------------------------
Polygon shclip(Polygon argp, int a, int b, int c, int d){ // Use Sutherland-Hodgman Algorithm to clip the polygon.
    // Clip edges of each polygon
    vector<Point> vecbuff;
    vector<Point> vecp = argp.getparr();
    Point buffp1, buffp2;
    Line buffl;
    int x1;
    int y1;
    int x2;
    int y2;
    // Clipping edge: x=a
    if (vecp.size() >= 1) { // Due to the for loop setting, we need to use a different method if vecp.size<1.
        if (a <= vecp[0].getx()) { // S-h Algorithm preprarition: If first point is inside, add it.
            vecbuff.push_back(argp.getp(0));
        }
        for (int i = 0; i < vecp.size()-1; i++) {
            buffp1 = vecp[i];
            buffp2 = vecp[i+1];
            x1 = buffp1.getx();
            x2 = buffp2.getx();
            if (a<=x1 && a<=x2) { // S-h Algorithm case 1
                vecbuff.push_back(buffp2);
            }else if (a<=x1 && x2<a) { // S-h Algorithm case 2
                double z1 = buffp1.getz();
                double z2 = buffp2.getz();
                double z_tmp;
                buffl.setp(buffp1, buffp2);
                int ya = buffl.cal_inter_wlr(a);
                if (abs(buffl.getslope()) <= 1) {
                    z_tmp = z1-(z1-z2)*(x1-a)/(x1-x2);
                }else {
                    y1 = buffp1.gety();
                    y2 = buffp2.gety();
                    z_tmp = z1-(z1-z2)*(y1-ya)/(y1-y2);
                }
                buffp2.set(a, ya);
                buffp2.zbuffer(z_tmp);
                if (buffp1 != buffp2) { // This part is to avoid pushing consecutive same vertices into one polygon.
                    vecbuff.push_back(buffp2);
                }
            }else if (x1<a && x2>=a) { // S-h Algorithm case 4
                double z1 = buffp1.getz();
                double z2 = buffp2.getz();
                double z_tmp;
                buffl.setp(buffp1, buffp2);
                int ya = buffl.cal_inter_wlr(a);
                if (abs(buffl.getslope()) <= 1) {
                    z_tmp = z1-(z1-z2)*(x1-a)/(x1-x2);
                }else {
                    y1 = buffp1.gety();
                    y2 = buffp2.gety();
                    z_tmp = z1-(z1-z2)*(y1-ya)/(y1-y2);
                }
                buffp1.set(a, ya);
                buffp1.zbuffer(z_tmp);
                if (buffp1 != buffp2) { // Avoiding pushing consecutive same vertices into one polygon
                    vecbuff.push_back(buffp1);
                }
                vecbuff.push_back(buffp2);
            } // else: S-h Algorithm case 3, just ignore
        }
        vecp.clear();
        vecp = vecbuff;
        vecbuff.clear();
        if (vecp.back() != vecp[0]) {
            vecp.push_back(vecp[0]);
        }
    }else vecp.clear();
    // Clipping edge: y=b
    if (vecp.size() >= 1) { // Due to the for loop setting, we need to use a different method if vecp.size<1.
        if (b <= vecp[0].gety()) { // S-h Algorithm preprarition: If first point is inside, add it.
            vecbuff.push_back(vecp[0]);
        }
        for (int i = 0; i < vecp.size()-1; i++) {
            buffp1 = vecp[i];
            buffp2 = vecp[i+1];
            y1 = buffp1.gety();
            y2 = buffp2.gety();
            if (b<=y1 && b<=y2) { // S-h Algorithm case 1
                vecbuff.push_back(buffp2);
            }else if (b<=y1 && y2<b) { // S-h Algorithm case 2
                double z1 = buffp1.getz();
                double z2 = buffp2.getz();
                double z_tmp;
                buffl.setp(buffp1, buffp2);
                int xb = buffl.cal_inter_wtb(b);
                if (abs(buffl.getslope()) >= 1) {
                    z_tmp = z1-(z1-z2)*(y1-b)/(y1-y2);
                }else {
                    x1 = buffp1.getx();
                    x2 = buffp2.getx();
                    z_tmp = z1-(z1-z2)*(x1-xb)/(x1-x2);
                }
                buffp2.set(xb, b);
                buffp2.zbuffer(z_tmp);
                if (buffp1 != buffp2) {
                    vecbuff.push_back(buffp2);
                }
            }else if (y1<b && b<=y2) { // S-h Algorithm case 4
                double z1 = buffp1.getz();
                double z2 = buffp2.getz();
                double z_tmp;
                buffl.setp(buffp1, buffp2);
                int xb = buffl.cal_inter_wtb(b);
                if (abs(buffl.getslope()) >= 1) {
                    z_tmp = z1-(z1-z2)*(y1-b)/(y1-y2);
                }else {
                    x1 = buffp1.getx();
                    x2 = buffp2.getx();
                    z_tmp = z1-(z1-z2)*(x1-xb)/(x1-x2);
                }
                buffp1.set(xb, b);
                buffp1.zbuffer(z_tmp);
                if (buffp1 != buffp2) {
                    vecbuff.push_back(buffp1);
                }
                vecbuff.push_back(buffp2);
            } // else: S-h Algorithm case 3, just ignore
        }
        vecp.clear();
        vecp = vecbuff;
        vecbuff.clear();
        if (vecp.back() != vecp[0]) {
            vecp.push_back(vecp[0]);
        }
    } else vecp.clear();
    // Clipping edge: x=c
    if (vecp.size() >= 1) { // Due to the for loop setting, we need to use a different method if vecp.size<1.
        if (vecp[0].getx() <= c) { // S-h Algorithm preprarition: If first point is inside, add it.
            vecbuff.push_back(vecp[0]);
        }
        for (int i = 0; i < vecp.size()-1; i++) {
            buffp1 = vecp[i];
            buffp2 = vecp[i+1];
            x1 = buffp1.getx();
            x2 = buffp2.getx();
            if (x1<=c && x2<=c) { // S-h Algorithm case 1
                vecbuff.push_back(buffp2);
            }else if (x1<=c && c<x2) { // S-h Algorithm case 2
                double z1 = buffp1.getz();
                double z2 = buffp2.getz();
                double z_tmp;
                buffl.setp(buffp1, buffp2);
                int yc = buffl.cal_inter_wlr(c);
                if (abs(buffl.getslope()) <= 1) {
                    z_tmp = z1-(z1-z2)*(x1-c)/(x1-x2);
                }else {
                    y1 = buffp1.gety();
                    y2 = buffp2.gety();
                    z_tmp = z1-(z1-z2)*(y1-yc)/(y1-y2);
                }
                buffp2.set(c, yc);
                buffp2.zbuffer(z_tmp);
                if (buffp1 != buffp2) {
                    vecbuff.push_back(buffp2);
                }
            }else if (c<x1 && x2<=c) { // S-h Algorithm case 4
                double z1 = buffp1.getz();
                double z2 = buffp2.getz();
                double z_tmp;
                buffl.setp(buffp1, buffp2);
                int yc = buffl.cal_inter_wlr(c);
                if (abs(buffl.getslope()) <= 1) {
                    z_tmp = z1-(z1-z2)*(x1-c)/(x1-x2);
                }else {
                    y1 = buffp1.gety();
                    y2 = buffp2.gety();
                    z_tmp = z1-(z1-z2)*(y1-yc)/(y1-y2);
                }
                buffp1.set(c, yc);
                buffp1.zbuffer(z_tmp);
                if (buffp1 != buffp2) {
                    vecbuff.push_back(buffp1);
                }
                vecbuff.push_back(buffp2);
            } // else: S-h Algorithm case 3, just ignore
        }
        vecp.clear();
        vecp = vecbuff;
        vecbuff.clear();
        if (vecp.back() != vecp[0]) {
            vecp.push_back(vecp[0]);
        }
    }else vecp.clear();
    // Clipping edge: y=d
    if (vecp.size() >= 1) {
        if (vecp[0].gety() <= d) { // S-h Algorithm preprarition: If first point is inside, add it.
            vecbuff.push_back(vecp[0]);
        }
        for (int i = 0; i < vecp.size()-1; i++) {
            buffp1 = vecp[i];
            buffp2 = vecp[i+1];
            y1 = buffp1.gety();
            y2 = buffp2.gety();
            if (y1<=d && y2<=d) { // S-h Algorithm case 1
                vecbuff.push_back(buffp2);
            }else if (y1<=d && d<y2) { // S-h Algorithm case 2
                double z1 = buffp1.getz();
                double z2 = buffp2.getz();
                double z_tmp;
                buffl.setp(buffp1, buffp2);
                int xd = buffl.cal_inter_wtb(d);
                if (abs(buffl.getslope()) >= 1) {
                    z_tmp = z1-(z1-z2)*(y1-d)/(y1-y2);
                }else {
                    x1 = buffp1.getx();
                    x2 = buffp2.getx();
                    z_tmp = z1-(z1-z2)*(x1-xd)/(x1-x2);
                }
                buffp2.set(xd, d);
                buffp2.zbuffer(z_tmp);
                if (buffp1 != buffp2) {
                    vecbuff.push_back(buffp2);
                }
            }else if (d<y1 && y2<=d) { // S-h Algorithm case 4
                double z1 = buffp1.getz();
                double z2 = buffp2.getz();
                double z_tmp;
                buffl.setp(buffp1, buffp2);
                int xd = buffl.cal_inter_wtb(d);
                if (abs(buffl.getslope()) >= 1) {
                    z_tmp = z1-(z1-z2)*(y1-d)/(y1-y2);
                }else {
                    x1 = buffp1.getx();
                    x2 = buffp2.getx();
                    z_tmp = z1-(z1-z2)*(x1-xd)/(x1-x2);
                }
                buffp1.set(xd, d);
                buffp1.zbuffer(z_tmp);
                if (buffp1 != buffp2) {
                    vecbuff.push_back(buffp1);
                }
                vecbuff.push_back(buffp2);
            } // else: S-h Algorithm case 3, just ignore
        }
        vecp.clear();
        vecp = vecbuff;
        vecbuff.clear();
        if (vecp.back() != vecp[0]) {
            vecp.push_back(vecp[0]);
        }
    }else vecp.clear();
    argp.setp(vecp);
    vecp.clear();
    return argp;
}

//---------------------------------------------------------------------------------------------------
void readfile(string input, vector<Vertex> *vecver, vector<Face> *vecf){
    // Start to read file.
    ifstream infile(input.c_str());
    if (!infile) {
        cout<<"Cannot open your file, please check your input path."<<endl;
        abort();
    }
    bool storev = false;
    bool storef = false;
    string str;
    vector<double> buff;
    while (infile) {
        infile>>str;
        if (str.compare("v") == 0) {
            storev = true;
        }else if (str.compare("f") == 0){
            storef = true;
        }else if (storev == true){ // Add vertex to its vector.
            buff.push_back(atof(str.c_str()));
            if (buff.size() == 3) {
                Vertex ver;
                ver.set(buff[0], buff[1], buff[2]);
                vecver->push_back(ver);
                storev = false;
                buff.clear();
            }
        }else if (storef == true){ // Add face to its vector.
            buff.push_back(atoi(str.c_str()));
            if (buff.size() == 3) {
                Face fac;
                fac.p1 = (int)buff[0] - 1;
                fac.p2 = (int)buff[1] - 1;
                fac.p3 = (int)buff[2] - 1;
                vecf->push_back(fac);
                storef = false;
                buff.clear();
            }
        }
    }
    infile.close();
}

//---------------------------------------------------------------------------------------------------
void outfile(string output){
    ofstream out(output.c_str());
    string line = "";
    if (!out) {
        cout<<"Cannot write an output file, please check your output path."<<endl;
    }
    out<<setheader()<<endl;
    int height = 501;
    int width = 501;
    cout<<setheader()<<endl;
    for (int i = height - 1; i >= 0; i--) {
        for (int j = 0; j< width; j++) {
            line += xpmpixel[i][j].color;
        }
        line = "\"" + line + "\"";
        if (i != 0) {
            line = line + ",";
        }
        out<<line<<endl;
        cout<<line<<endl;
        line.clear();
    }
    out<<setend()<<endl;
    cout<<setend()<<endl;
    out.close();
    
}

//---------------------------------------------------------------------------------------------------
string setheader() {
    stringstream ss;
    int w = 501;
    int h = 501;
    ss<<w;
    string intw;
    ss>>intw;
    ss.clear();
    ss<<h;
    string inth;
    ss>>inth;
    ss.clear();
    int colors = 1;
    string str = "";
    // Calculate 20 values for one kind of color.
    double step = (double)0xff / 20.0;
    vector<int> color_value;
    for (int i = 1 ; i <= 20; i++) {
        color_value.push_back(rnd(step*i));
    }
    if (!f[0].empty()) { // Set headers for red colors.
        colors += 20;
        string red = "";
        for (int i = 0; i < 20; i++) {
            char color_char = (char)(60+i);
            ss<<color_char;
            string color;
            ss>>color;
            ss.clear();
            red += "\"" + color + " c #";
            ss<<color_value[i];
            ss>>hex>>color;
            ss.clear();
            if (color.length() == 1) {
                color = "0" + color;
            }
            red += color + "0000\",\n";
        }
        str += red;
    }
    if (!f[1].empty()) { // Set headers for green colors.
        colors += 20;
        string green = "";
        for (int i = 0; i < 20; i++) {
            char color_char = (char)(80+i);
            string color;
            ss<<color_char;
            ss>>color;
            ss.clear();
            green += "\"" + color + " c #00";
            ss<<color_value[i];
            ss>>hex>>color;
            ss.clear();
            if (color.length() == 1) {
                color = "0" + color;
            }
            green += color + "00\",\n";
        }
        str += green;
    }
    if (!f[2].empty()) { // Set headers for blue colors.
        colors += 20;
        string blue = "";
        for (int i = 0; i < 20; i++) {
            char color_char = (char)(100+i);
            string color;
            ss<<color_char;
            ss>>color;
            ss.clear();
            blue += "\"" + color + " c #0000";
            ss<<color_value[i];
            ss>>hex>>color;
            ss.clear();
            if (color.length() == 1) {
                color = "0" + color;
            }
            blue += color + "\",\n";
        }
        str += blue;
    }
    ss<<dec<<colors;
    string color;
    ss>>dec>>color;
    str = "/* XPM */\nstatic char *CG_hw5[] = {\n/* width height num_colors chars_per_pixel */\n\"" + intw + " " + inth + " " + color + " 1\",\n/* colors */\n\"- c #000000\",\n" + str;
    str += "/* pixels */";
    return str;
}

//---------------------------------------------------------------------------------------------------
string setend() {
    string str = "};";
    return str;
}

//---------------------------------------------------------------------------------------------------
int rnd(double arg){
    if (arg >= 0) {
        return (int)(arg + 0.5);
    }else{
        return (int)(arg - 0.5);
    }
}

//---------------------------------------------------------------------------------------------------
void optana(int argc, char * const argv[]){
    // analyze input option and set default options
    int opt;
    while ((opt = getopt(argc, argv, "j:k:o:p:x:y:z:X:Y:Z:q:r:w:Q:R:W:u:v:U:V:F:B:Pf:g:i:h"))!= -1) {
        switch (opt) {
            case 'j':
                j = atoi(optarg);
                break;
            case 'k':
                k = atoi(optarg);
                break;
            case 'o':
                o = atoi(optarg);
                break;
            case 'p':
                p = atoi(optarg);
                break;
            case 'x':
                x = atof(optarg);
                break;
            case 'y':
                y = atof(optarg);
                break;
            case 'z':
                z = atof(optarg);
                break;
            case 'X':
                X  = atof(optarg);
                break;
            case 'Y':
                Y = atof(optarg);
                break;
            case 'Z':
                Z = atof(optarg);
                break;
            case 'q':
                q = atof(optarg);
                break;
            case 'r':
                r = atof(optarg);
                break;
            case 'w':
                w = atof(optarg);
                break;
            case 'Q':
                Q = atof(optarg);
                break;
            case 'R':
                R = atof(optarg);
                break;
            case 'W':
                W = atof(optarg);
                break;
            case 'u':
                u = atof(optarg);
                break;
            case 'v':
                v = atof(optarg);
                break;
            case 'U':
                U = atof(optarg);
                break;
            case 'V':
                V = atof(optarg);
                break;
            case 'F':
                F = atof(optarg);
                break;
            case 'B':
                B = atof(optarg);
                break;
            case 'f':
                f[0] = optarg;
                break;
            case 'g':
                f[1] = optarg;
                break;
            case 'i':
                f[2] = optarg;
                break;
            case 'P':
                P = true;
                break;
            case 'h':
                help();
                abort();
                break;
            default:cout<<"Your input commands may not be correct. Please enter -h for help."<<endl;
                abort();
                break;
        }
    }
    if (j > o) {
        int tmp;
        tmp = j;
        j = o;
        o = tmp;
    }
    if (k > p) {
        int tmp;
        tmp = k;
        k = p;
        p = tmp;
    }
    if ((o>500) || (p>500)) {
        cout<<"The viewport cannot be larger than xpm image."<<endl;
        abort();
    }
    if (f[0].empty() && f[1].empty() && f[2].empty()) {
        f[0] = "./bound-sprellpsd.smf";
    }
    if (F < 0) {
        cout<<"Value of F must be positive."<<endl;
        abort();
    }
    if (B > 0) {
        cout<<"Value of B must be negative."<<endl;
        abort();
    }
}

//------------------------------------------------------------------------------------------
void help(){
    cout<<"This program implements z-buffering with near and far plane clipping. The output XPM image has up to 60 colors, which are determined by distance from screen. Deeper the vertex inside the screen, darker the pixel displays."<<endl;
    cout<<"Output XPM image will be saved in the same working path with name of out.xpm"<<endl;
    cout<<"Following are the features and input options available of this program."<<endl;
    cout<<"[-f] Followed by the name of the first SMF model. Its base color is \"ff0000\" (red)."<<endl;
    cout<<"[-g] Followed by the name of the second SMF model. Its base color is \"00ff00\" (green)."<<endl;
    cout<<"[-i] Followed by the name of the third SMF model. Its base color is \"0000ff\" (blue)."<<endl;
    cout<<"[-j] The next argument is an integer lower bound in the x dimension of the view port in screen coordinates"<<endl;
    cout<<"[-k] The next argument is an integer in lower bound in the y dimension of the view port in screen coordinates"<<endl;
    cout<<"The next argument is an integer in upper bound in the x dimension of the view port in screen coordinates"<<endl;
    cout<<"[-p] The next argument is an integer in upper bound in the y dimension of the view port in screen coordinates"<<endl;
    cout<<"[-x] floating point x of Projection Reference Point (PRP) in VRC coordinates"<<endl;
    cout<<"[-y] floating point y of Projection Reference Point (PRP) in VRC coordinates"<<endl;
    cout<<"[-z] floating point z of Projection Reference Point (PRP) in VRC coordinates"<<endl;
    cout<<"[-X] floating point x of View Reference Point (VRP) in world coordinates"<<endl;
    cout<<"[-Y] floating point y of View Reference Point (VRP) in world coordinates"<<endl;
    cout<<"[-Z] floating point z of View Reference Point (VRP) in world coordinates"<<endl;
    cout<<"[-q] floating point x of View Plane Normal vector (VPN) in world coordinates"<<endl;
    cout<<"[-r] floating point y of View Plane Normal vector (VPN) in world coordinates"<<endl;
    cout<<"[-w] floating point z of View Plane Normal vector (VPN) in world coordinates"<<endl;
    cout<<"[-Q] floating point x of View Up Vector (VUP) in world coordinates"<<endl;
    cout<<"[-R] floating point y of View Up Vector (VUP) in world coordinates"<<endl;
    cout<<"[-W] floating point z of View UP Vector (VUP) in world coordinates"<<endl;
    cout<<"[-u] floating point u min of the VRC window in VRC coordinates"<<endl;
    cout<<"[-v] floating point v min of the VRC window in VRC coordinates"<<endl;
    cout<<"[-U] floating point u max of the VRC window in VRC coordinates"<<endl;
    cout<<"[-V] floating point v max of the VRC window in VRC coordinates"<<endl;
    cout<<"[-P] Use parallel projection. If this flag is not present, use perspective projection."<<endl;
    cout<<"[-F] Followed by the floating point coordinate of the Front (Near) plane in VRC coordinates."<<endl;
    cout<<"[-B] Followed by the floating point coordinate of the Back (Far) plane in VRC coordinates."<<endl;
    cout<<"The default values are: "<<endl;
    cout<<"./CG_hw5 -f bound-sprellpsd.smf -j 0 -k 0 -o 500 -p 500 -x 0.0 -y 0.0 -z 1.0 -X 0.0 -Y 0.0 -Z 0.0 -q 0.0 -r 0.0 -w -1.0 -Q 0.0 -R 1.0 -W 0.0 -u -0.7 -v -0.7 -U 0.7 -V 0.7 -F 0.6 -B -0.6 > out.xpm"<<endl;
}

//------------------------------------------------------------------------------------------
bool scancmp(scanline scan_in1, scanline scan_in2) {
    // This function is for sorting of scanline elements.
    return scan_in1.x<scan_in2.x;
}