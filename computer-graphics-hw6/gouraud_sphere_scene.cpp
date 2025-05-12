//
//  sphere_scene.c
//  Rasterizer
//
//
#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <cmath>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int     gNumVertices1 = 0;    // Number of 3D vertices.
int     gNumTriangles1 = 0;    // Number of triangles.
int* gIndexBuffer1 = NULL; // Vertex indices for the triangles.
unsigned char image1[512][512][3];
float z_buffer1[512][512];

struct Vec3 {
	float x, y, z;
};

Vec3 operator-(const Vec3& v) {
	return { -v.x, -v.y, -v.z };
}

Vec3* gVertexBuffer1 = nullptr;
Vec3* gColorBuffer = nullptr;

Vec3* gNormalBuffer = nullptr;
Vec3 operator+(const Vec3& a, const Vec3& b) { return { a.x + b.x, a.y + b.y, a.z + b.z }; }
Vec3 operator-(const Vec3& a, const Vec3& b) { return { a.x - b.x, a.y - b.y, a.z - b.z }; }
Vec3 operator*(const Vec3& a, float s) { return { a.x * s, a.y * s, a.z * s }; }

Vec3 normalize1(const Vec3& v) {
	float len = sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
	return { v.x / len, v.y / len, v.z / len };
}

Vec3 cross1(const Vec3& a, const Vec3& b) {
	return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
}

void compute_vertex_normals1() {
	gNormalBuffer = new Vec3[gNumVertices1];
	for (int i = 0; i < gNumVertices1; ++i)
		gNormalBuffer[i] = { 0, 0, 0 };

	for (int i = 0; i < gNumTriangles1; ++i) {
		int k0 = gIndexBuffer1[3 * i + 0];
		int k1 = gIndexBuffer1[3 * i + 1];
		int k2 = gIndexBuffer1[3 * i + 2];

		Vec3 v0 = gVertexBuffer1[k0];
		Vec3 v1 = gVertexBuffer1[k1];
		Vec3 v2 = gVertexBuffer1[k2];

		Vec3 n = normalize1(cross1(v2 - v0, v1 - v0));
		gNormalBuffer[k0] = gNormalBuffer[k0] + n;
		gNormalBuffer[k1] = gNormalBuffer[k1] + n;
		gNormalBuffer[k2] = gNormalBuffer[k2] + n;
	}

	for (int i = 0; i < gNumVertices1; ++i)
		gNormalBuffer[i] = normalize1(gNormalBuffer[i]);
}

float dot1(const Vec3& a, const Vec3& b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 compute_lighting1(const Vec3& pos, const Vec3& normal) {
	Vec3 ka = { 0.0f, 1.0f, 0.0f };
	Vec3 kd = { 0.0f, 0.5f, 0.0f };
	Vec3 ks = { 0.5f, 0.5f, 0.5f };
	float p = 32.0f;
	float Ia = 0.2f;
	Vec3 light_pos = { -4, 4, -3 };
	Vec3 light_color = { 1, 1, 1 };

	Vec3 n = normalize1(normal);
	Vec3 l = normalize1(light_pos - pos);
	Vec3 v = normalize1(-pos); // assuming eye at origin
	Vec3 h = normalize1(l + v);

	float diff = std::max(0.0f, dot1(n, l));
	float spec = powf(std::max(0.0f, dot1(n, h)), p);

	Vec3 color = ka * Ia + kd * diff + ks * spec;
	color.x = powf(std::min(1.0f, color.x), 1.0f / 2.2f);
	color.y = powf(std::min(1.0f, color.y), 1.0f / 2.2f);
	color.z = powf(std::min(1.0f, color.z), 1.0f / 2.2f);
	return color;
}

Vec3 model_transform1(const Vec3& v) {
	return { v.x * 2, v.y * 2, v.z * 2 - 7 };
}

void compute_vertex_colors1() {
	gColorBuffer = new Vec3[gNumVertices1];
	for (int i = 0; i < gNumVertices1; ++i) {
		Vec3 world_pos = model_transform1(gVertexBuffer1[i]);
		gColorBuffer[i] = compute_lighting1(world_pos, gNormalBuffer[i]);
	}
}

Vec3 perspective_transform1(const Vec3& v) {
	float l = -0.1f, r = 0.1f, b = -0.1f, t = 0.1f, n = -0.1f, f = -1000.0f;

	float x = (2 * n * v.x) / (r - l);
	float y = (2 * n * v.y) / (t - b);
	float z = (f + n + 2 * f * n / v.z) / (f - n);
	return { x, y, z };
}

Vec3 viewport_transform1(const Vec3& v) {
	int nx = 512, ny = 512;
	return {
		(v.x + 1.0f) * 0.5f * nx,
		(1.0f - v.y) * 0.5f * ny,
		v.z
	};
}

Vec3 transform_vertex1(const Vec3& v) {
	Vec3 model = model_transform1(v);
	Vec3 clip = perspective_transform1(model);
	Vec3 ndc = { clip.x / -model.z, clip.y / -model.z, clip.z };
	Vec3 screen = viewport_transform1(ndc);
	return screen;
}

void clear_buffers1() {
	for (int y = 0; y < 512; ++y) {
		for (int x = 0; x < 512; ++x) {
			image1[y][x][0] = image1[y][x][1] = image1[y][x][2] = 0;
			z_buffer1[y][x] = 1e9f;
		}
	}
}

void draw_triangle1(Vec3 v0, Vec3 v1, Vec3 v2, Vec3 c0, Vec3 c1, Vec3 c2) {
	// Bounding box
	float min_x = std::min(std::min(v0.x, v1.x), v2.x);
	float max_x = std::max(std::max(v0.x, v1.x), v2.x);
	float min_y = std::min(std::min(v0.y, v1.y), v2.y);
	float max_y = std::max(std::max(v0.y, v1.y), v2.y);

	int minX = std::max(0, (int)std::floor(min_x));
	int maxX = std::min(511, (int)std::ceil(max_x));
	int minY = std::max(0, (int)std::floor(min_y));
	int maxY = std::min(511, (int)std::ceil(max_y));

	float denom = (v1.y - v2.y) * (v0.x - v2.x) + (v2.x - v1.x) * (v0.y - v2.y);

	for (int y = minY; y <= maxY; ++y) {
		for (int x = minX; x <= maxX; ++x) {
			float px = x + 0.5f;
			float py = y + 0.5f;

			float alpha = ((v1.y - v2.y) * (px - v2.x) + (v2.x - v1.x) * (py - v2.y)) / denom;
			float beta = ((v2.y - v0.y) * (px - v2.x) + (v0.x - v2.x) * (py - v2.y)) / denom;
			float gamma = 1.0f - alpha - beta;

			if (alpha >= 0 && beta >= 0 && gamma >= 0) {
				float z = alpha * v0.z + beta * v1.z + gamma * v2.z;
				if (z < z_buffer1[y][x]) {
					z_buffer1[y][x] = z;
					Vec3 c = c0 * alpha + c1 * beta + c2 * gamma;
					image1[y][x][0] = (unsigned char)(std::min(1.0f, c.x) * 255);
					image1[y][x][1] = (unsigned char)(std::min(1.0f, c.y) * 255);
					image1[y][x][2] = (unsigned char)(std::min(1.0f, c.z) * 255);
				}
			}
		}
	}
}

void render_scene1() {
	for (int i = 0; i < gNumTriangles1; ++i) {
		int k0 = gIndexBuffer1[3 * i + 0];
		int k1 = gIndexBuffer1[3 * i + 1];
		int k2 = gIndexBuffer1[3 * i + 2];

		Vec3 v0 = transform_vertex1(gVertexBuffer1[k0]);
		Vec3 v1 = transform_vertex1(gVertexBuffer1[k1]);
		Vec3 v2 = transform_vertex1(gVertexBuffer1[k2]);

		Vec3 c0 = gColorBuffer[k0];
		Vec3 c1 = gColorBuffer[k1];
		Vec3 c2 = gColorBuffer[k2];

		draw_triangle1(v0, v1, v2, c0, c1, c2);
	}
}

void save_image_bmp1(const char* filename) {
	FILE* f;
	int w = 512, h = 512;
	int filesize = 54 + 3 * w * h;
	unsigned char bmpfileheader[14] = {
		'B','M', filesize & 255, (filesize >> 8) & 255, (filesize >> 16) & 255, (filesize >> 24) & 255,
		0,0, 0,0, 54,0,0,0
	};
	unsigned char bmpinfoheader[40] = {
		40,0,0,0, w & 255, (w >> 8) & 255, (w >> 16) & 255, (w >> 24) & 255,
		h & 255, (h >> 8) & 255, (h >> 16) & 255, (h >> 24) & 255,
		1,0, 24,0, 0,0,0,0,
		0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0
	};

	f = fopen(filename, "wb");
	fwrite(bmpfileheader, 1, 14, f);
	fwrite(bmpinfoheader, 1, 40, f);

	for (int y = h - 1; y >= 0; y--) {
		for (int x = 0; x < w; x++) {
			unsigned char r = image1[y][x][0];
			unsigned char g = image1[y][x][1];
			unsigned char b = image1[y][x][2];
			unsigned char color[3] = { b, g, r };
			fwrite(color, 1, 3, f);
		}
	}
	fclose(f);
}

void create_scene1()
{
	int width = 32;
	int height = 16;

	float theta, phi;
	int t;

	gNumVertices1 = (height - 2) * width + 2;
	gNumTriangles1 = (height - 2) * (width - 1) * 2;
	gVertexBuffer1 = new Vec3[gNumVertices1];

	// TODO: Allocate an array for gNumVertices vertices.

	gIndexBuffer1 = new int[3 * gNumTriangles1];

	t = 0;
	for (int j = 1; j < height - 1; ++j)
	{
		for (int i = 0; i < width; ++i)
		{
			theta = (float)j / (height - 1) * M_PI;
			phi = (float)i / (width - 1) * M_PI * 2;

			float   x = sinf(theta) * cosf(phi);
			float   y = cosf(theta);
			float   z = -sinf(theta) * sinf(phi);

			// TODO: Set vertex t in the vertex array to {x, y, z}.

			gVertexBuffer1[t++] = { x, y, z };
		}
	}

	// TODO: Set vertex t in the vertex array to {0, 1, 0}.

	gVertexBuffer1[t++] = { 0.0f, 1.0f, 0.0f };

	// TODO: Set vertex t in the vertex array to {0, -1, 0}.

	gVertexBuffer1[t++] = { 0.0f, -1.0f, 0.0f };

	t = 0;
	for (int j = 0; j < height - 3; ++j)
	{
		for (int i = 0; i < width - 1; ++i)
		{
			gIndexBuffer1[t++] = j * width + i;
			gIndexBuffer1[t++] = (j + 1) * width + (i + 1);
			gIndexBuffer1[t++] = j * width + (i + 1);
			gIndexBuffer1[t++] = j * width + i;
			gIndexBuffer1[t++] = (j + 1) * width + i;
			gIndexBuffer1[t++] = (j + 1) * width + (i + 1);
		}
	}
	for (int i = 0; i < width - 1; ++i)
	{
		gIndexBuffer1[t++] = (height - 2) * width;
		gIndexBuffer1[t++] = i;
		gIndexBuffer1[t++] = i + 1;
		gIndexBuffer1[t++] = (height - 2) * width + 1;
		gIndexBuffer1[t++] = (height - 3) * width + (i + 1);
		gIndexBuffer1[t++] = (height - 3) * width + i;
	}

	// The index buffer has now been generated. Here's how to use to determine
	// the vertices of a triangle. Suppose you want to determine the vertices
	// of triangle i, with 0 <= i < gNumTriangles. Define:
	//
	// k0 = gIndexBuffer[3*i + 0]
	// k1 = gIndexBuffer[3*i + 1]
	// k2 = gIndexBuffer[3*i + 2]
	//
	// Now, the vertices of triangle i are at positions k0, k1, and k2 (in that
	// order) in the vertex array (which you should allocate yourself at line
	// 27).
	//
	// Note that this assumes 0-based indexing of arrays (as used in C/C++,
	// Java, etc.) If your language uses 1-based indexing, you will have to
	// add 1 to k0, k1, and k2.
}

int main2() {
	create_scene1();
	compute_vertex_normals1();
	compute_vertex_colors1();
	clear_buffers1();
	render_scene1();
	save_image_bmp1("gouraud .bmp");
	return 0;
}
