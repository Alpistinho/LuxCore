/***************************************************************************
 * Copyright 1998-2018 by authors (see AUTHORS.txt)                        *
 *                                                                         *
 *   This file is part of LuxCoreRender.                                   *
 *                                                                         *
 * Licensed under the Apache License, Version 2.0 (the "License");         *
 * you may not use this file except in compliance with the License.        *
 * You may obtain a copy of the License at                                 *
 *                                                                         *
 *     http://www.apache.org/licenses/LICENSE-2.0                          *
 *                                                                         *
 * Unless required by applicable law or agreed to in writing, software     *
 * distributed under the License is distributed on an "AS IS" BASIS,       *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 * See the License for the specific language governing permissions and     *
 * limitations under the License.                                          *
 ***************************************************************************/

#include <map>
#include <vector>
#include <string>

#include <boost/format.hpp>

#include "luxrays/core/exttrianglemesh.h"
#include "slg/shapes/simplify.h"
#include "slg/scene/scene.h"
#include "slg/utils/harlequincolors.h"

using namespace std;
using namespace luxrays;
using namespace slg;

//------------------------------------------------------------------------------
//
// The following code is based on Sven Forstmann's quadric mesh simplification
// code (https://github.com/sp4cerat/Fast-Quadric-Mesh-Simplification)
// and heavily modified for LuxCoreRender

/////////////////////////////////////////////
//
// Mesh Simplification Tutorial
//
// (C) by Sven Forstmann in 2014
//
// License : MIT
// http://opensource.org/licenses/MIT
//
//https://github.com/sp4cerat/Fast-Quadric-Mesh-Simplification
//
// 5/2016: Chris Rorden created minimal version for OSX/Linux/Windows compile

class SymetricMatrix {
public:
	// Constructor
	SymetricMatrix(const float c = 0) {
		for (u_int i = 0; i < 10; ++i)
			m[i] = c;
	}

	SymetricMatrix(
			const float m11, const float m12, const float m13, const float m14,
			const float m22, const float m23, const float m24,
			const float m33, const float m34,
			const float m44) {
		m[0] = m11;
		m[1] = m12;
		m[2] = m13;
		m[3] = m14;
		m[4] = m22;
		m[5] = m23;
		m[6] = m24;
		m[7] = m33;
		m[8] = m34;
		m[9] = m44;
	}

	// Make plane
	SymetricMatrix(const float a, const float b, const float c, const float d) {
		m[0] = a * a;
		m[1] = a * b;
		m[2] = a * c;
		m[3] = a * d;
		m[4] = b * b;
		m[5] = b * c;
		m[6] = b * d;
		m[7 ] = c * c;
		m[8 ] = c * d;
		m[9 ] = d * d;
	}

	float operator[](int c) const {
		return m[c];
	}

	// Determinant
	float det(
			const u_int a11, const u_int a12, const u_int a13,
			const u_int a21, const u_int a22, const u_int a23,
			const u_int a31, const u_int a32, const u_int a33) const {
		const float det = m[a11] * m[a22] * m[a33] + m[a13] * m[a21] * m[a32] + m[a12] * m[a23] * m[a31]
				- m[a13] * m[a22] * m[a31] - m[a11] * m[a23] * m[a32] - m[a12] * m[a21] * m[a33];
		return det;
	}

	const SymetricMatrix operator+(const SymetricMatrix &n) const {
		return SymetricMatrix(
				m[0] + n[0], m[1] + n[1], m[2] + n[2], m[3] + n[3],
				m[4] + n[4], m[5] + n[5], m[6] + n[6],
				m[7] + n[7], m[8] + n[8],
				m[9] + n[9]);
	}

	SymetricMatrix& operator+=(const SymetricMatrix& n) {
		m[0] += n[0];
		m[1] += n[1];
		m[2] += n[2];
		m[3] += n[3];
		m[4] += n[4];
		m[5] += n[5];
		m[6] += n[6];
		m[7] += n[7];
		m[8] += n[8];
		m[9] += n[9];

		return *this;
	}

	float m[10];
};


class Simplify {
private:
	struct Triangle {
		u_int v[3];
		Normal geometryN;
		float err[4];
		bool deleted, dirty;
	};

	struct Vertex {
		Point p;
		u_int tstart, tcount;
		SymetricMatrix q;
		bool border;
	};

	struct Ref {
		int tid, tvertex;
	};
	
	vector<Triangle> triangles;
	vector<Vertex> vertices;
	vector<Ref> refs;

public:
	Simplify(const ExtTriangleMesh &srcMesh) {
		const u_int vertCount = srcMesh.GetTotalVertexCount();
		const u_int triCount = srcMesh.GetTotalTriangleCount();
		const Point *verts = srcMesh.GetVertices();
		const luxrays::Triangle *tris = srcMesh.GetTriangles();

		vertices.resize(vertCount);
		for (u_int i = 0; i < vertCount; ++i)
			vertices[i].p = verts[i];

		triangles.resize(triCount);
		for (u_int i = 0; i < triCount; ++i) {
			triangles[i].v[0] = tris[i].v[0];
			triangles[i].v[1] = tris[i].v[1];
			triangles[i].v[2] = tris[i].v[2];
		}
	}

	~Simplify() {
	}
	
	ExtTriangleMesh *GetExtMesh() const {
		const u_int vertCount = vertices.size();
		const u_int triCount = triangles.size();
		Point *newVertices = ExtTriangleMesh::AllocVerticesBuffer(vertCount);
		luxrays::Triangle *newTris = ExtTriangleMesh::AllocTrianglesBuffer(triCount);
		
		for (u_int i = 0; i < vertCount; ++i)
			newVertices[i] = vertices[i].p;
		
		for (u_int i = 0; i < triCount; ++i) {
			newTris[i].v[0] = triangles[i].v[0];
			newTris[i].v[1] = triangles[i].v[1];
			newTris[i].v[2] = triangles[i].v[2];
		}
		
		return new ExtTriangleMesh(vertCount, triCount, newVertices, newTris);
	}
	
	void Decimate(const float errorScale) {
		// Init
		for (u_int i = 0; i < triangles.size(); ++i)
			triangles[i].deleted = false;

		// Main iteration loop
		u_int deletedTriangles = 0;
		vector<bool> deleted0, deleted1;
		for (u_int iteration = 0;; ++iteration) {
			if (triangles.size() <= 1)
				break;

			// Update mesh constantly
			UpdateMesh(iteration);

			// Clear dirty flag
			for (u_int i = 0; i < triangles.size(); ++i)
				triangles[i].dirty = false;

			//
			// All triangles with edges below the threshold will be removed
			//
			// The following numbers works well for most models.
			// If it does not, try to adjust the 3 parameters
			//
			const float threshold = errorScale * .001f;
			
			SDL_LOG("Simplify iteration: " << iteration);

			// Remove vertices & mark deleted triangles
			for (u_int i = 0; i < triangles.size(); ++i) {
				Triangle &t = triangles[i];

				if (t.err[3] > threshold)
					continue;
				if (t.deleted)
					continue;
				if (t.dirty)
					continue;

				for (u_int j = 0; j < 3; ++j) {
					if (t.err[j] < threshold) {
						u_int i0 = t.v[j];
						Vertex &v0 = vertices[i0];
						
						const u_int i1 = t.v[(j + 1) % 3];
						Vertex &v1 = vertices[i1];

						// Border check
						if (v0.border != v1.border)
							continue;

						// Compute vertex to collapse to
						Point p;
						CalculateError(i0, i1, p);

						deleted0.resize(v0.tcount); // normals temporarily
						deleted1.resize(v1.tcount); // normals temporarily

						// don't remove if flipped
						if (Flipped(p, i0, i1, v0, v1, deleted0))
							continue;
						if (Flipped(p, i1, i0, v1, v0, deleted1))
							continue;

						// Not flipped, so remove edge
						v0.p = p;
						v0.q = v1.q + v0.q;
						u_int tstart = refs.size();

						UpdateTriangles(i0, v0, deleted0, deletedTriangles);
						UpdateTriangles(i0, v1, deleted1, deletedTriangles);

						const u_int tcount = refs.size() - tstart;

						if (tcount <= v0.tcount) {
							// save ram
							if (tcount)memcpy(&refs[v0.tstart], &refs[tstart], tcount * sizeof (Ref));
						} else
							// append
							v0.tstart = tstart;

						v0.tcount = tcount;
						break;
					}
				}
			}
			if (deletedTriangles <= 0)break;
			deletedTriangles = 0;
		}

		// Clean up mesh
		CompactMesh();
	}

private:
	// Check if a triangle flips when this edge is removed
	bool Flipped(const Point &p, const u_int i0, const u_int i1,
			const Vertex &v0, const Vertex &v1,
			vector<bool> &deleted) {
		for (u_int k = 0; k < v0.tcount; ++k) {
			Triangle &t = triangles[refs[v0.tstart + k].tid];

			if (t.deleted)
				continue;

			const u_int s = refs[v0.tstart + k].tvertex;
			const u_int id1 = t.v[(s + 1) % 3];
			const u_int id2 = t.v[(s + 2) % 3];

			// Delete ?
			if (id1 == i1 || id2 == i1) {
				deleted[k] = true;
				continue;
			}

			const Vector d1 = Normalize(vertices[id1].p - p);
			const Vector d2 = Normalize(vertices[id2].p - p);

			if (AbsDot(d1, d2) > .999f)
				return true;

			const Normal geometryN(Normalize(Cross(d1, d2)));
			deleted[k] = false;
			if (Dot(geometryN, t.geometryN) < .2f)
				return true;
		}

		return false;
	}


	// Update triangle connections and edge error after a edge is collapsed
	void UpdateTriangles(const u_int i0, const Vertex &v, const  vector<bool> &deleted,
			u_int &deletedTriangles) {
		Point p;

		for (u_int k = 0; k < v.tcount; ++k) {
			Ref &r = refs[v.tstart + k];
			Triangle &t = triangles[r.tid];

			if (t.deleted)
				continue;

			if (deleted[k]) {
				t.deleted = 1;
				deletedTriangles++;
				continue;
			}

			t.v[r.tvertex] = i0;
			t.dirty = true;
			t.err[0] = CalculateError(t.v[0], t.v[1], p);
			t.err[1] = CalculateError(t.v[1], t.v[2], p);
			t.err[2] = CalculateError(t.v[2], t.v[0], p);
			t.err[3] = Min(t.err[0], Min(t.err[1], t.err[2]));

			refs.push_back(r);
		}
	}

	// Compact triangles, compute edge error and build reference list
	void UpdateMesh(const u_int iteration) {
		if (iteration > 0) {
			// compact triangles
			int dst = 0;
			for (u_int i = 0; i < triangles.size(); ++i)
				if (!triangles[i].deleted)
					triangles[dst++] = triangles[i];

			triangles.resize(dst);
		}

		//
		// Init Quadrics by Plane & Edge Errors
		//
		// required at the beginning ( iteration == 0 )
		// recomputing during the simplification is not required,
		// but mostly improves the result for closed meshes
		//
		if (iteration == 0) {
			for (u_int i = 0; i < vertices.size(); ++i)
				vertices[i].q = SymetricMatrix(0.0);

			for (u_int i = 0; i < triangles.size(); ++i) {
				Triangle &t = triangles[i];

				Point p[3];
				p[0] = vertices[t.v[0]].p;
				p[1] = vertices[t.v[1]].p;
				p[2] = vertices[t.v[2]].p;
				
				const Normal geometryN(Normalize(Cross(p[1] - p[0], p[2] - p[0])));
				t.geometryN = geometryN;

				const SymetricMatrix sm(geometryN.x, geometryN.y, geometryN.z,
						-Dot(Vector(geometryN), Vector(p[0])));
				vertices[t.v[0]].q = vertices[t.v[0]].q + sm;
				vertices[t.v[1]].q = vertices[t.v[1]].q + sm;
				vertices[t.v[2]].q = vertices[t.v[2]].q + sm;
			}

			for (u_int i = 0; i < triangles.size(); ++i) {
				// Calc Edge Error
				Triangle &t = triangles[i];

				Point p;
				t.err[0] = CalculateError(t.v[0], t.v[1], p);
				t.err[1] = CalculateError(t.v[1], t.v[2], p);
				t.err[2] = CalculateError(t.v[2], t.v[0], p);
				
				t.err[3] = Min(t.err[0], Min(t.err[1], t.err[2]));
			}
		}

		// Init Reference ID list
		for (u_int i = 0; i < vertices.size(); ++i) {
			vertices[i].tstart = 0;
			vertices[i].tcount = 0;
		}

		for (u_int i = 0; i < triangles.size(); ++i) {
			Triangle &t = triangles[i];

			vertices[t.v[0]].tcount++;
			vertices[t.v[1]].tcount++;
			vertices[t.v[2]].tcount++;
		}

		u_int tstart = 0;
		for (u_int i = 0; i < vertices.size(); ++i) {
			Vertex &v = vertices[i];

			v.tstart = tstart;
			tstart += v.tcount;
			v.tcount = 0;
		}

		// Write References
		refs.resize(triangles.size() * 3);
		for (u_int i = 0; i < triangles.size(); ++i) {
			Triangle &t = triangles[i];

			for (u_int j = 0; j < 3; ++j) {
				Vertex &v = vertices[t.v[j]];

				refs[v.tstart + v.tcount].tid = i;
				refs[v.tstart + v.tcount].tvertex = j;

				v.tcount++;
			}
		}

		// Identify boundary : vertices[].border=0,1
		if (iteration == 0) {
			for (u_int i = 0; i < vertices.size(); ++i)
				vertices[i].border = false;

			vector<u_int> vcount, vids;
			for (u_int i = 0; i < vertices.size(); ++i) {
				Vertex &v = vertices[i];
				vcount.clear();
				vids.clear();

				for (u_int j = 0; j < v.tcount; ++j) {
					int k = refs[v.tstart + j].tid;
					Triangle &t = triangles[k];

					for (u_int k = 0; k < 3; ++k) {
						u_int ofs = 0;
						u_int id = t.v[k];

						while (ofs < vcount.size()) {
							if (vids[ofs] == id)
								break;

							ofs++;
						}

						if (ofs == vcount.size()) {
							vcount.push_back(1);
							vids.push_back(id);
						} else
							vcount[ofs]++;
					}
				}
				
				for (u_int j = 0; j < vcount.size(); ++j) {
					if (vcount[j] == 1)
						vertices[vids[j]].border = 1;
				}
			}
		}
	}

	// Finally compact mesh before exiting
	void CompactMesh() {
		u_int dst = 0;

		for (u_int i = 0; i < vertices.size(); ++i)
			vertices[i].tcount = 0;

		for (u_int i = 0; i < triangles.size(); ++i) {
			if (!triangles[i].deleted) {
				Triangle &t = triangles[i];
				triangles[dst++] = t;

				vertices[t.v[0]].tcount = 1;
				vertices[t.v[1]].tcount = 1;
				vertices[t.v[2]].tcount = 1;
			}
		}
		triangles.resize(dst);

		dst = 0;
		for (u_int i = 0; i < vertices.size(); ++i) {
			if (vertices[i].tcount) {
				vertices[i].tstart = dst;
				vertices[dst].p = vertices[i].p;
				dst++;
			}
		}

		for (u_int i = 0; i < triangles.size(); ++i) {
			Triangle &t = triangles[i];

			t.v[0] = vertices[t.v[0]].tstart;
			t.v[1] = vertices[t.v[1]].tstart;
			t.v[2] = vertices[t.v[2]].tstart;
		}
		vertices.resize(dst);
	}

	// Error between vertex and Quadric
	float VertexError(const SymetricMatrix &q, const float x, const float y, const float z) {
		return q[0] * x * x + 2 * q[1] * x * y + 2 * q[2] * x * z + 2 * q[3] * x + q[4] * y * y
				+ 2 * q[5] * y * z + 2 * q[6] * y + q[7] * z * z + 2 * q[8] * z + q[9];
	}

	// Error for one edge
	float CalculateError(const u_int id_v1, const u_int id_v2, Point &pResult) {
		// Compute interpolated vertex
		const SymetricMatrix q = vertices[id_v1].q + vertices[id_v2].q;
		const bool border = vertices[id_v1].border && vertices[id_v2].border;

		float error = 0.f;
		const float det = q.det(0, 1, 2, 1, 4, 5, 2, 5, 7);
		if ((det != 0.f) && !border) {
			// q_delta is invertible
			pResult.x = -1 / det * (q.det(1, 2, 3, 4, 5, 6, 5, 7, 8)); // vx = A41/det(q_delta)
			pResult.y = 1 / det * (q.det(0, 2, 3, 1, 5, 6, 2, 7, 8)); // vy = A42/det(q_delta)
			pResult.z = -1 / det * (q.det(0, 1, 3, 1, 4, 6, 2, 5, 8)); // vz = A43/det(q_delta)

			error = VertexError(q, pResult.x, pResult.y, pResult.z);
		} else {
			// det = 0 -> try to find best result
			const Point &p1 = vertices[id_v1].p;
			const Point &p2 = vertices[id_v2].p;
			const Point p3 = (p1 + p2) / 2;

			const float error1 = VertexError(q, p1.x, p1.y, p1.z);
			const float error2 = VertexError(q, p2.x, p2.y, p2.z);
			const float error3 = VertexError(q, p3.x, p3.y, p3.z);
			error = Min(error1, Min(error2, error3));

			if (error1 == error)
				pResult = p1;
			if (error2 == error)
				pResult = p2;
			if (error3 == error)
				pResult = p3;
		}

		return error;
	}
};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

SimplifyShape::SimplifyShape(luxrays::ExtTriangleMesh *srcMesh, const float errorScale) {
	SDL_LOG("Simplify shape " << srcMesh->GetName() << " with error scale: " << errorScale);

	const float startTime = WallClockTime();

	Simplify simplify(*srcMesh);
	simplify.Decimate(errorScale);
	mesh = simplify.GetExtMesh();

	SDL_LOG("Subdivided shape from " << srcMesh->GetTotalTriangleCount() << " to " << mesh->GetTotalTriangleCount() << " faces");

	// For some debugging
	//mesh->Save("debug.ply");
	
	const float endTime = WallClockTime();
	SDL_LOG("Simplify time: " << (boost::format("%.3f") % (endTime - startTime)) << "secs");
}

SimplifyShape::~SimplifyShape() {
	if (!refined)
		delete mesh;
}

ExtTriangleMesh *SimplifyShape::RefineImpl(const Scene *scene) {
	return mesh;
}
