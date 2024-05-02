//Cells raycast
{
	PVSRayHit tmpHit;
	float cellMin[3];
	float cellMax[3];
	
	size_t cellPos[3] = {startPos[0], startPos[1], startPos[2]};
	while(cellPos[0] != hitPos[0] && cellPos[1] != hitPos[1] && cellPos[2] != hitPos[2]) {
		
		PVSCell* cell = db->cells + (cellPos[0] + cellPos[1] * db->cellsX + cellPos[2] * db->cellsX * db->cellsY);

		cell->visMesh[rayHit.hitMeshId / 8] |= 1 << (rayHit.hitMeshId % 8);
		
		//X plane test
		if(ray.dir[0] < 0 && cellPos[0] > 0) {
			cellPos[0] -= 1;
			_getCellBounds(db, cellPos, cellMin, cellMax);
			
			if(raycastCube(&tmpHit, ray, cellMin, cellMax)) continue;
			cellPos[0] += 1;
		} else if(ray.dir > 0 && cellPos[0] < db->cellsX - 1) {
			cellPos[0] += 1;
			_getCellBounds(db, cellPos, cellMin, cellMax);
			
			if(raycastCube(&tmpHit, ray, cellMin, cellMax)) continue;
			cellPos[0] -= 1;
		}
		
		//Y plane test
		if(ray.dir[1] < 0 && cellPos[1] > 0) {
			cellPos[1] -= 1;
			_getCellBounds(db, cellPos, cellMin, cellMax);
			
			if(raycastCube(&tmpHit, ray, cellMin, cellMax)) continue;
			cellPos[1] += 1;
		} else if(ray.dir[1] > 0 && cellPos[1] < db->cellsY - 1) {
			cellPos[1] += 1;
			_getCellBounds(db, cellPos, cellMin, cellMax);
			
			if(raycastCube(&tmpHit, ray, cellMin, cellMax)) continue;
			cellPos[1] -= 1;
		}
		
		//Z plane test
		if(ray.dir[2] < 0 && cellPos[2] > 0) {
			cellPos[2] -= 1;
			_getCellBounds(db, cellPos, cellMin, cellMax);
			
			if(raycastCube(&tmpHit, ray, cellMin, cellMax)) continue;
			cellPos[2] += 1;
		} else if(ray.dir[2] > 0 && cellPos[2] < db->cellsY - 1) {
			cellPos[2] += 1;
			_getCellBounds(db, cellPos, cellMin, cellMax);
			
			if(raycastCube(&tmpHit, ray, cellMin, cellMax)) continue;
			cellPos[2] -= 1;
		}
		
		//not hits?
		break;
	}
}
void _pvsProcessRaycastResult(PVSdb* db, const int startMeshId, const PVSRayHit rayHit, const PVSRay ray) {
	size_t startPos[3];
	_getCellPos(db, ray.start, startPos);
	
	//Used for hit position detection
	float hitDistance = rayHit.hitDistance;
	if(!rayHit.hit) hitDistance = sqrtf(pow2f(db->max[0] - db->min[0]) + pow2f(db->max[1] - db->min[1]) + pow2f(db->max[2] - db->min[2]));
	
	size_t hitPos[3];
	float tmp[3] = {
		ray.start[0] + ray.dir[0] * hitDistance,
		ray.start[1] + ray.dir[1] * hitDistance,
		ray.start[2] + ray.dir[2] * hitDistance
	};
	_getCellPos(db, tmp, hitPos);
	
	PVSRayHit tmpHit;
	size_t minP[3] = {mins(startPos[0], hitPos[0]), mins(startPos[1], hitPos[1]), mins(startPos[2], hitPos[2])};
	size_t maxP[3] = {maxs(startPos[0], hitPos[0]), maxs(startPos[1], hitPos[1]), maxs(startPos[2], hitPos[2])};
	
	for(size_t x=minP[0]; x<=maxP[0]; x++) {
		for(size_t y=minP[1]; y<=maxP[1]; y++) {

			for(size_t z=minP[2]; z<=maxP[2]; z++) {
				//Ray vs cell collision
				size_t cellPos[3] = {x, y, z};

				float cellMin[3];
				float cellMax[3];
				_getCellBounds(db, cellPos, cellMin, cellMax);

				tmpHit.hit = false;
				raycastCube(&tmpHit, ray, cellMin, cellMax);

				if(tmpHit.hit) {
					PVSCell* cell = db->cells + (x + y * db->cellsX + z * db->cellsX * db->cellsY);

					if(rayHit.hit && !rayHit.hitBackFace) cell->visMesh[rayHit.hitMeshId / 8] |= 1 << (rayHit.hitMeshId % 8);
					cell->visMesh[startMeshId / 8] |= 1 << (startMeshId % 8);
				}
			}
		}
	}
}

size_t pvsCompute(PVSdb* db, PVSModel mdl, float raysDensity, float adaptRaysCountDistance) {
	srand(time(NULL));
	float volume = (db->max[0] - db->min[0]) * (db->max[1] - db->min[1]) * (db->max[2] - db->min[2]);
	size_t raysCount = (int) roundf(volume * raysDensity);
	
	for(size_t i=0; i<raysCount; i++) {
		
		PVSRay ray;
		PVSRayHit rayHit;
		rayHit.hit = false;
		
		//Starting position
		for(int axis=0; axis<3; axis++)
			ray.start[axis] = (db->max[axis] - db->min[axis]) * randf() + db->min[axis];
		
		//Starting direction
		float rotX = randf() * M_PI  - M_PI / 2; //from -90 to +90
		float rotY = randf() * M_PI * 2; //from 0 to 360
		
		ray.dir[0] = (float) sin(rotY) * cos(rotX);
		ray.dir[1] = (float) sin(rotX);
		ray.dir[2] = (float) cos(rotY) * cos(rotX);
		
		for(int t=0; t<cvector_size(mdl.meshes); t++) {
			PVSMesh mesh = mdl.meshes[t];
			
			if(rayAABBtest(ray, mesh.min, mesh.max)) {
				rayCast(&rayHit, mesh, ray);
			}
		}
		
		if(rayHit.hit && !rayHit.hitBackFace) _pvsProcessRaycastResult(db, rayHit, ray);
		
		if((i % (raysCount / 10)) == 0) {
			printf("%d %%\n", (int) (i * 100 / raysCount));
		}
	}
	
	return raysCount;
}

size_t _compareVis(char* vis1, char* vis2, size_t size) {
	size_t error = 0;
	
	for(size_t i=0; i<size; i++) {
		int cError = (vis1[i] ^ vis2[i]);// & (~vis2[i]);
		
		for(int b=0; b<8; b++) {
			error += (cError >> b) & 1;
		}
	}
	
	return error;
}

float _tryExtend(
	PVSdb* db, char* visMeshes, size_t visMeshesSize,
	bool* usedCells,
	int* errorArray,
	size_t x, size_t y, size_t z, 
	//float* errorArray,
	size_t* cellSize, 
	int axis
	) {
		
	assert(axis >= 0 && axis <= 2);
	
	if(axis == 0 && (x + cellSize[0]) >= db->cellsX) return INT_MAX;
	else if(axis == 1 && (y + cellSize[1]) >= db->cellsY) return INT_MAX;
	else if(axis == 2 && (z + cellSize[2]) >= db->cellsZ) return INT_MAX;
	
	size_t min[3] = {x, y, z};
	size_t max[3] = {x + cellSize[0], y + cellSize[1], z + cellSize[2]};
	
	size_t oldMin[3] = {min[0], min[1], min[2]};
	size_t oldMax[3] = {max[0], max[1], max[2]};
	
	if(axis == 0) {
		min[0] = max[0];
		max[0]++;
	} else if(axis == 1) {
		min[1] = max[1];
		max[1]++;
	} else if(axis == 2) {
		min[2] = max[2];
		max[2]++;
	}
	
	//memset(errorArray, 0, sizeof(float) * visMeshesSize * 8);
	int error = 0;
	
	for(size_t x2=min[0]; x2<max[0]; x2++) {
		for(size_t y2=min[1]; y2<max[1]; y2++) {
			for(size_t z2=min[2]; z2<max[2]; z2++) {
				
				if(usedCells[x2 + y2 * db->cellsX + z2 * db->cellsX * db->cellsY]) return INT_MAX;
				
				PVSCell* cell = db->cellGrid + (x2 + y2 * db->cellsX + z2 * db->cellsX * db->cellsY);
				
				error += _compareVis(visMeshes, cell->visMesh/*, errorArray*/, visMeshesSize);
				
				//if(memcmp(visMeshes, cell->visMesh, visMeshesSize)) return false;
				/*for(size_t x3=oldMin[0]; x3<oldMax[0]; x3++) {
					for(size_t y3=oldMin[1]; y3<oldMax[1]; y3++) {
						for(size_t z3=oldMin[2]; z3<oldMax[2]; z3++) {
							
							PVSCell* cell2 = db->cellGrid + (x3 + y3 * db->cellsX + z3 * db->cellsX * db->cellsY);
							
							for(size_t i=0; i<visMeshesSize*8; i++) {
								int vis1 = (cell->visMesh[i / 8] >> (i%8)) & 1;
								int vis2 = (cell2->visMesh[i / 8] >> (i%8)) & 1;
								int vis3 = (visMeshes[i / 8] >> (i%8)) & 1;
								
								error += (vis1 == 1 && vis2 == 0 && vis3 == 0) ? errorArray[i] : 0;
							}
							
						}
					}
				}*/
				
			}
		}
	}
	
	return error;
}

cvector_vector_type(PVSCell) pvsJoinCells(PVSdb* db, int maxError) {
	
	cvector_vector_type(PVSCell) cells = NULL;
	
	size_t visMeshesSize = (db->meshes / 8) + ((db->meshes % 8) > 0 ? 1 : 0);
	char* tmpVisMeshes = malloc(visMeshesSize);
	assert(tmpVisMeshes);
	
	int* errorArray = malloc(visMeshesSize * 8 * sizeof(int));
	assert(errorArray);
	
	bool* usedCells = calloc(db->cellsX * db->cellsY * db->cellsZ, sizeof(bool));
	assert(usedCells);
	
	for(size_t x=0; x<db->cellsX; x++) {
		for(size_t y=0; y<db->cellsY; y++) {
			for(size_t z=0; z<db->cellsZ; z++) {
				
				//Cell is already used, continue search
				if(usedCells[x + y * db->cellsX + z * db->cellsX * db->cellsY]) continue;
				
				//Found unused cell, try to extend it
				PVSCell* cell = db->cellGrid + (x + y * db->cellsX + z * db->cellsX * db->cellsY);
				memcpy(tmpVisMeshes, cell->visMesh, visMeshesSize);
				
				//mark it as used
				usedCells[x + y * db->cellsX + z * db->cellsX * db->cellsY] = true;
				
				size_t cellSize[3] = {1, 1, 1};
				int currentError = 0;
				
				while(true) {
					//size_t cellSizeOld[3] = {cellSize[0], cellSize[1], cellSize[2]};
					memset(errorArray, 0, sizeof(int) * visMeshesSize * 8);
					
					for(size_t x2=0; x2<cellSize[0]; x2++) {
						for(size_t y2=0; y2<cellSize[1]; y2++) {
							for(size_t z2=0; z2<cellSize[2]; z2++) {
								
								PVSCell* cell2 = db->cellGrid + ((x + x2) + (y + y2) * db->cellsX + (z + z2) * db->cellsX * db->cellsY);
								
								for(size_t i=0; i<visMeshesSize*8; i++) 
									errorArray[i] += ((cell2->visMesh[i / 8]) & (1 << (i%8))) == 0;
							}
						}
					}*/
					
					int err1 = _tryExtend(db, tmpVisMeshes, visMeshesSize, usedCells, errorArray, x, y, z, cellSize, 0);
					int err2 = _tryExtend(db, tmpVisMeshes, visMeshesSize, usedCells, errorArray, x, y, z, cellSize, 1);
					int err3 = _tryExtend(db, tmpVisMeshes, visMeshesSize, usedCells, errorArray, x, y, z, cellSize, 2);
					
					/*float err1f = (float)err1 / cellSize[1] / cellSize[2];
					float err2f = (float)err2 / cellSize[0] / cellSize[2];
					float err3f = (float)err3 / cellSize[1] / cellSize[0];*/
					
					//Select axis with minimal error
					if(err1 != INT_MAX && (currentError + err1f) <= maxError && err1 <= err2 && err1 <= err3) {
						cellSize[0]++;
						//currentError += err1;
					} if(err2 != INT_MAX && (currentError + err2f) <= maxError && err2 <= err1 && err2 <= err3) {
						cellSize[1]++;
						//currentError += err2;
					} if(err3 != INT_MAX && (currentError + err3f) <= maxError && err3 <= err1 && err3 <= err2) {
						cellSize[2]++;
						//currentError += err3;
					} else break;
					
					//if(cellSizeOld[0] == cellSize[0] && cellSizeOld[1] == cellSize[1] && cellSizeOld[2] == cellSize[2]) break;
					
					//mark as used, update visibility data
					for(size_t x2=0; x2<cellSize[0]; x2++) {
						for(size_t y2=0; y2<cellSize[1]; y2++) {
							for(size_t z2=0; z2<cellSize[2]; z2++) {
								
								usedCells[(x + x2) + (y + y2) * db->cellsX + (z + z2) * db->cellsX * db->cellsY] = true;
								
								PVSCell* cell2 = db->cellGrid + ((x + x2) + (y + y2) * db->cellsX + (z + z2) * db->cellsX * db->cellsY);
								
								for(size_t i=0; i<visMeshesSize; i++) {
									tmpVisMeshes[i] |= cell2->visMesh[i];
								}
							}
						}
					}
				}
				
				/*for(size_t x2=0; x2<cellSize[0]; x2++) {
					for(size_t y2=0; y2<cellSize[1]; y2++) {
						for(size_t z2=0; z2<cellSize[2]; z2++) {
							
							usedCells[(x + x2) + (y + y2) * db->cellsX + (z + z2) * db->cellsX * db->cellsY] = true;
							
							PVSCell* cell2 = db->cellGrid + ((x + x2) + (y + y2) * db->cellsX + (z + z2) * db->cellsX * db->cellsY);
							
							for(size_t i=0; i<visMeshesSize; i++) {
								tmpVisMeshes[i] |= cell2->visMesh[i];
							}
						}
					}
				}*/
				
				PVSCell newCell = {0};
				
				newCell.min[0] = x;
				newCell.min[1] = y;
				newCell.min[2] = z;
				newCell.max[0] = x + cellSize[0];
				newCell.max[1] = y + cellSize[1];
				newCell.max[2] = z + cellSize[2];
				
				char* newVisMesh = malloc(visMeshesSize);
				assert(newVisMesh);
				memcpy(newVisMesh, tmpVisMeshes, visMeshesSize);
				newCell.visMesh = newVisMesh;
				
				cvector_push_back(cells, newCell);
				
				//Reset search
				x = 0;
				y = 0;
				z = 0;
				
				//printf("Cells: %d / %d\n", (int)cvector_size(cells), (int)(db->cellsX * db->cellsY * db->cellsZ));
			}
		}
	}
	
	return cells;
}