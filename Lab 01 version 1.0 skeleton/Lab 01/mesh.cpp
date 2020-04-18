#include "mesh.h"

#ifdef _WIN32
#include <Windows.h>
#include "GL\glut.h"
#define M_PI 3.141592654
#elif __APPLE__
#include <OpenGL/gl.h>
#include <GLUT/GLUT.h>
#endif

#include "math.h"
#include <string>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "mesh.h"
#include <map>
#include <queue>
#include <iomanip>
#include <set>
#include <queue>   
using namespace std;



void myObjType::draw() {

	glEnable(GL_LIGHTING);

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

	glPushMatrix();
	double longestSide = 0.0;
	for (int i = 0; i < 3; i++)
		if ((lmax[i] - lmin[i]) > longestSide)
			longestSide = (lmax[i] - lmin[i]);
	glScalef(4.0 / longestSide, 4.0 / longestSide, 4.0 / longestSide);
	glTranslated(-(lmin[0] + lmax[0]) / 2.0, -(lmin[1] + lmax[1]) / 2.0, -(lmin[2] + lmax[2]) / 2.0);
	for (int i = 1; i <= tcount; i++)
	{
		glBegin(GL_POLYGON);
// uncomment the following after you computed the normals
		glNormal3dv(nlist[i]);
		for (int j = 0; j < 3; j++)
			glVertex3dv(vlist[tlist[i][j]]);
		glEnd();
	
	}
	glDisable(GL_LIGHTING);

	glPopMatrix();
}

void myObjType::writeFile(char* filename)
{
	std::ofstream outfile;
	outfile.open(filename, std::ios::out);
	std::string strbuff;
	for (int i = 1; i <= vcount; ++i) {
		outfile << "v " << vlist[i][0] << " " << vlist[i][1] << " " << vlist[i][2] << std::endl;
	}
	outfile << std::endl;
	for (int i = 1; i <= tcount; ++i) {
		outfile << "f " << tlist[i][0] << " " << tlist[i][1] << " " << tlist[i][2] << std::endl;
	}
	outfile.close();
}

void myObjType::readFile(char* filename)
{
	cout << "Opening " << filename << endl;
	ifstream inFile;
	inFile.open(filename);
	if (!inFile.is_open()) {
		cout << "We cannot find your file " << filename << endl;
		exit(1);
	}

	string line;
	int i, j;
	bool firstVertex = 1;
	double currCood;

	while (getline(inFile, line))
	{
		if ((line[0] == 'v' || line[0] == 'f') && line[1] == ' ')
		{
			if (line[0] == 'v')
			{
				vcount++;
				i = 1;
				const char* linec = line.data();
				for (int k = 0; k < 3; k++) { // k is 0,1,2 for x,y,z
					while (linec[i] == ' ') i++;
					j = i;
					while (linec[j] != ' ') j++;
					currCood = vlist[vcount][k] = atof(line.substr(i, j - i).c_str());
					if (firstVertex) 
						lmin[k] = lmax[k] = currCood;
					else {
						if (lmin[k] > currCood)
							lmin[k] = currCood;
						if (lmax[k] < currCood)
							lmax[k] = currCood;
					}
					i = j;
				}

				firstVertex = 0;
			}
			if (line[0] == 'f')
			{
				tcount++;
				i = 1;
				const char* linec = line.data();
				for (int k = 0; k < 3; k++) {
					while (linec[i] == ' ') i++;
					j = i;
					while (linec[j] != ' ' && linec[j] != '\\') j++;
					tlist[tcount][k] = atof(line.substr(i, j - i).c_str());
					i = j;
					fnlist[tcount][k] = 0;
					while (linec[j] != ' ') j++;

				}

			}


		}
	}

	// We suggest you to compute the normals here
	computeTrianglesNormals();
	findFNext();
    cout << "No. of vertices: " << vcount << endl;
    cout << "No. of triangles: " << tcount << endl;
    computeStat();
}

int myObjType::org(OrTri ot) {
	int v = ver(ot);
	int id = idx(ot);
	if (v < 3) {
		return tlist[id][v];
	}
	else {
		return tlist[id][(v+1)%3];
	}
}

int myObjType::dest(OrTri ot) {
	int v = ver(ot);
	int id = idx(ot);
	if (v < 3) {
		return tlist[id][(v+1)%3];
	}
	else {
		return tlist[id][v % 3];
	}

}

void myObjType::findFNext(){
	for (int i = 1; i <= tcount; ++i) {
		for (int j = 0; j < 3; ++j) {
			int oriTri = makeOrTri(i, j);
			int start = org(oriTri), end = dest(oriTri);
			std::pair <int, int>key1(start, end);	
			std::map<pair <int, int>, int>::iterator it1 = hashMap.find(key1);
			if (it1 != hashMap.end()) {
				int found_oriTri = hashMap.at(key1);
				int id = idx(found_oriTri);
				int v = ver(found_oriTri);
				fnlist[i][j]= found_oriTri;
				fnlist[id][v%3] = oriTri;
				hashMap.erase(it1);
			}
			else{
				std::pair <int, int>key2(end, start);
				std::map<pair <int, int>, int>::iterator it2 = hashMap.find(key2);
				if (it2 != hashMap.end()) {
					int found_oriTri = hashMap.at(key2);
					int id = idx(found_oriTri);
					int v = ver(found_oriTri);
					fnlist[i][j] = found_oriTri;
					fnlist[id][v % 3] = oriTri;
					hashMap.erase(it2);
				}
				else {
					hashMap.emplace(key1, oriTri);
				}
			}
		}
	}

	for (int i = 1; i <= tcount; ++i) {
		for (int j = 0; j < 3; ++j) {
			if (fnlist[i][j] == 0) {
				fnlist[i][j] = makeOrTri(i, j);
			}
		}
	}
	
}

void myObjType::computeTrianglesNormals() {
	for (int i = 1; i <= tcount; ++i) {
		double ax = vlist[tlist[i][1]][0] - vlist[tlist[i][0]][0];
		double ay = vlist[tlist[i][1]][1] - vlist[tlist[i][0]][1];
		double az = vlist[tlist[i][1]][2] - vlist[tlist[i][0]][2];

		double bx = vlist[tlist[i][2]][0] - vlist[tlist[i][0]][0];
		double by = vlist[tlist[i][2]][1] - vlist[tlist[i][0]][1];
		double bz = vlist[tlist[i][2]][2] - vlist[tlist[i][0]][2];

		double nx = (ay * bz) - (az * by);
		double ny = (az * bx) - (ax * bz);
		double nz = (ax * by) - (ay * bx);

		double length = sqrt(nx * nx + ny * ny + nz * nz);
		if (length > 0.0) {
			nx /= length;
			ny /= length;
			nz /= length;
		}
		nlist[i][0] = nx;
		nlist[i][1] = ny;
		nlist[i][2] = nz;
	}
}



void myObjType::computeAngles() {
	for (int i = 1; i <= tcount; ++i) {
		for (int j = 0; j < 3; ++j) {
			int first;
			int second;
			switch (j) {
				case 0: {
					first = 1;
					second = 2;
					break;
				}
				case 1: {
					first = 0;
					second = 2;
					break;
				}
				case 2: {
					first = 1;
					second = 0;
					break;
				}

			}
			double ax = vlist[tlist[i][first]][0] - vlist[tlist[i][j]][0];
			double ay = vlist[tlist[i][first]][1] - vlist[tlist[i][j]][1];
			double az = vlist[tlist[i][first]][2] - vlist[tlist[i][j]][2];
			double alength = sqrt(ax * ax + ay * ay + az * az);

			double bx = vlist[tlist[i][second]][0] - vlist[tlist[i][j]][0];
			double by = vlist[tlist[i][second]][1] - vlist[tlist[i][j]][1];
			double bz = vlist[tlist[i][second]][2] - vlist[tlist[i][j]][2];
			double blength = sqrt(bx * bx + by * by + bz * bz);

			float alpha = acos((ax * bx + ay * by + az * bz) / (alength * blength)) * 180 / M_PI;
			int bucket = floor(alpha / 10);
			++statMaxAngle[bucket];
			++statMinAngle[bucket];
			if (alpha > maxAngle) {
				maxAngle = alpha;
			}
			if (alpha < minAngle) {
				minAngle = alpha;
			}
		}
	}
}

void myObjType::computeComponents() {
	std::set<int> triangles;
	for (int i = 1; i <= tcount; ++i) {
		triangles.insert(i);
	}
	components = 0;
	while (!triangles.empty()) {
		++components;
		std::priority_queue<int> queue;
		std::set<int> compontentTriangles;
		int first = *triangles.begin();
		queue.push(first);
		compontentTriangles.insert(first);
		while (!queue.empty()) {
			int index = queue.top();
			queue.pop();
			if (compontentTriangles.find(idx(fnlist[index][0]))==compontentTriangles.end()&&triangles.find(idx(fnlist[index][0]))!=triangles.end()) {
				queue.push(idx(fnlist[index][0]));
				compontentTriangles.insert(idx(fnlist[index][0]));
			}
			if (compontentTriangles.find(idx(fnlist[index][1]))==compontentTriangles.end()&&triangles.find(idx(fnlist[index][1]))!=triangles.end()) {
				queue.push(idx(fnlist[index][1]));
				compontentTriangles.insert(idx(fnlist[index][1]));
			}
			if (compontentTriangles.find(idx(fnlist[index][2]))==compontentTriangles.end()&&triangles.find(idx(fnlist[index][2]))!=triangles.end()) {
				queue.push(idx(fnlist[index][2]));
				compontentTriangles.insert(idx(fnlist[index][2]));
			}
			triangles.erase(index);
		}
	}
}

void myObjType::OriTriPrint() {
	std::ofstream outfile;
	outfile.open("log.txt", std::ios::out);
	std::string strbuff;
	for (int i = 1; i <= tcount; ++i) {
		outfile << "Triangle "<<i<<":  OriTri:"<<makeOrTri(i,0)<<" List:"<< fnlist[i][0] << " " << fnlist[i][1] << " " << fnlist[i][2] << std::endl;
	}
	outfile.close();
}

void myObjType::computeStat()
{
	int i;
	computeAngles();
	computeComponents();
    cout << "Min. angle = " << minAngle << endl;
    cout << "Max. angle = " << maxAngle << endl;
	
	cout << "Statistics for Maximum Angles" << endl;
	for (i = 0; i < 18; i++)
		cout << statMaxAngle[i] << " ";
	cout << endl;
	cout << "Statistics for Minimum Angles" << endl;
	for (i = 0; i < 18; i++)
		cout << statMinAngle[i] << " ";
	cout << endl;
	cout << "Components: " << components<<endl;
	OriTriPrint();
}


void myObjType::load3DS(char* p_filename) {

}