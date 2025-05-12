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

int     gNumVertices = 0;    // Number of 3D vertices.
int     gNumTriangles = 0;    // Number of triangles.
int* gIndexBuffer = NULL; // Vertex indices for the triangles.
unsigned char image[512][512][3];
float z_buffer[512][512];

struct Vec3 {
	float x, y, z;
};
Vec3* gVertexBuffer = nullptr;

Vec3 normalize(const Vec3& v) {
	float len = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
	return { v.x / len, v.y / len, v.z / len };
}

float dot(const Vec3& a, const Vec3& b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 compute_lighting(const Vec3& pos, const Vec3& normal) {
	Vec3 ka = { 0, 1, 0 }, kd = { 0, 0.5, 0 }, ks = { 0.5, 0.5, 0.5 };
	Vec3 light_pos = { -4, 4, -3 };
	Vec3 Ia = { 0.2f, 0.2f, 0.2f };
	Vec3 I = { 1.0f, 1.0f, 1.0f };
	int p = 32;

	Vec3 L = normalize({ light_pos.x - pos.x, light_pos.y - pos.y, light_pos.z - pos.z });
	Vec3 N = normalize(normal);
	Vec3 V = normalize({ -pos.x, -pos.y, -pos.z }); // viewer at origin
	Vec3 H = normalize({ L.x + V.x, L.y + V.y, L.z + V.z });

	float diff = std::max(0.0f, dot(N, L));
	float spec = std::pow(std::max(0.0f, dot(N, H)), p);

	Vec3 color = {
		Ia.x * ka.x + I.x * (kd.x * diff + ks.x * spec),
		Ia.y * ka.y + I.y * (kd.y * diff + ks.y * spec),
		Ia.z * ka.z + I.z * (kd.z * diff + ks.z * spec)
	};
	return color;
}

unsigned char gamma_correction(float c) {
	c = std::max(0.0f, std::min(1.0f, c));
	c = std::pow(c, 1.0f / 2.2f);
	return (unsigned char)(c * 255);
}

Vec3 model_transform(const Vec3& v) {
	return { v.x * 2, v.y * 2, v.z * 2 - 7 };
}

Vec3 perspective_transform(const Vec3& v) {
	float l = -0.1f, r = 0.1f, b = -0.1f, t = 0.1f, n = -0.1f, f = -1000.0f;

	float x = (2 * n * v.x) / (r - l);
	float y = (2 * n * v.y) / (t - b);
	float z = (f + n + 2 * f * n / v.z) / (f - n);
	return { x, y, z };
}

Vec3 viewport_transform(const Vec3& v) {
	int nx = 512, ny = 512;
	return {
		(v.x + 1.0f) * 0.5f * nx,
		(1.0f - v.y) * 0.5f * ny,
		v.z
	};
}

Vec3 transform_vertex(const Vec3& v) {
	Vec3 model = model_transform(v);
	Vec3 clip = perspective_transform(model);
	Vec3 ndc = { clip.x / -model.z, clip.y / -model.z, clip.z };
	Vec3 screen = viewport_transform(ndc);
	return screen;
}

void clear_buffers() {
	for (int y = 0; y < 512; ++y) {
		for (int x = 0; x < 512; ++x) {
			image[y][x][0] = image[y][x][1] = image[y][x][2] = 0;
			z_buffer[y][x] = 1e9f;
		}
	}
}

void draw_triangle(Vec3 v0, Vec3 v1, Vec3 v2, Vec3 world0, Vec3 world1, Vec3 world2) {
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

	Vec3 center = {
		(world0.x + world1.x + world2.x) / 3.0f,
		(world0.y + world1.y + world2.y) / 3.0f,
		(world0.z + world1.z + world2.z) / 3.0f
	};
	// ³ë¸Ö: (v1 - v0) x (v2 - v0)
	Vec3 a = {
		world2.x - world0.x,
		world2.y - world0.y,
		world2.z - world0.z
	};
	Vec3 b = {
		world1.x - world0.x,
		world1.y - world0.y,
		world1.z - world0.z
	};
	Vec3 normal = normalize({
		a.y * b.z - a.z * b.y,
		a.z * b.x - a.x * b.z,
		a.x * b.y - a.y * b.x
	});

	Vec3 color = compute_lighting(center, normal);
	unsigned char r = gamma_correction(color.x);
	unsigned char g = gamma_correction(color.y);
	unsigned char b_col = gamma_correction(color.z);

	for (int y = minY; y <= maxY; ++y) {
		for (int x = minX; x <= maxX; ++x) {
			float px = x + 0.5f;
			float py = y + 0.5f;

			float alpha = ((v1.y - v2.y) * (px - v2.x) + (v2.x - v1.x) * (py - v2.y)) / denom;
			float beta = ((v2.y - v0.y) * (px - v2.x) + (v0.x - v2.x) * (py - v2.y)) / denom;
			float gamma = 1.0f - alpha - beta;

			if (alpha >= 0 && beta >= 0 && gamma >= 0) {
				float z = alpha * v0.z + beta * v1.z + gamma * v2.z;
				if (z < z_buffer[y][x]) {
					z_buffer[y][x] = z;
					image[y][x][0] = r;
					image[y][x][1] = g;
					image[y][x][2] = b_col;
				}
			}
		}
	}
}

void render_scene() {
	for (int i = 0; i < gNumTriangles; ++i) {
		int k0 = gIndexBuffer[3 * i + 0];
		int k1 = gIndexBuffer[3 * i + 1];
		int k2 = gIndexBuffer[3 * i + 2];

		Vec3 world0 = model_transform(gVertexBuffer[k0]);
		Vec3 world1 = model_transform(gVertexBuffer[k1]);
		Vec3 world2 = model_transform(gVertexBuffer[k2]);

		Vec3 v0 = viewport_transform({ perspective_transform(world0).x / -world0.z, perspective_transform(world0).y / -world0.z, perspective_transform(world0).z });
		Vec3 v1 = viewport_transform({ perspective_transform(world1).x / -world1.z, perspective_transform(world1).y / -world1.z, perspective_transform(world1).z });
		Vec3 v2 = viewport_transform({ perspective_transform(world2).x / -world2.z, perspective_transform(world2).y / -world2.z, perspective_transform(world2).z });

		draw_triangle(v0, v1, v2, world0, world1, world2);
	}
}

void save_image_bmp(const char* filename) {
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
			unsigned char r = image[y][x][0];
			unsigned char g = image[y][x][1];
			unsigned char b = image[y][x][2];
			unsigned char color[3] = { b, g, r };
			fwrite(color, 1, 3, f);
		}
	}
	fclose(f);
}

void create_scene()
{
	int width = 32;
	int height = 16;

	float theta, phi;
	int t;

	gNumVertices = (height - 2) * width + 2;
	gNumTriangles = (height - 2) * (width - 1) * 2;
	gVertexBuffer = new Vec3[gNumVertices];

	// TODO: Allocate an array for gNumVertices vertices.

	gIndexBuffer = new int[3 * gNumTriangles];

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

			gVertexBuffer[t++] = { x, y, z };
		}
	}

	// TODO: Set vertex t in the vertex array to {0, 1, 0}.

	gVertexBuffer[t++] = { 0.0f, 1.0f, 0.0f };

	// TODO: Set vertex t in the vertex array to {0, -1, 0}.

	gVertexBuffer[t++] = { 0.0f, -1.0f, 0.0f };

	t = 0;
	for (int j = 0; j < height - 3; ++j)
	{
		for (int i = 0; i < width - 1; ++i)
		{
			gIndexBuffer[t++] = j * width + i;
			gIndexBuffer[t++] = (j + 1) * width + (i + 1);
			gIndexBuffer[t++] = j * width + (i + 1);
			gIndexBuffer[t++] = j * width + i;
			gIndexBuffer[t++] = (j + 1) * width + i;
			gIndexBuffer[t++] = (j + 1) * width + (i + 1);
		}
	}
	for (int i = 0; i < width - 1; ++i)
	{
		gIndexBuffer[t++] = (height - 2) * width;
		gIndexBuffer[t++] = i;
		gIndexBuffer[t++] = i + 1;
		gIndexBuffer[t++] = (height - 2) * width + 1;
		gIndexBuffer[t++] = (height - 3) * width + (i + 1);
		gIndexBuffer[t++] = (height - 3) * width + i;
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

int main1() {
	create_scene();
	clear_buffers();
	render_scene();
	save_image_bmp("flat.bmp");
	return 0;
}

