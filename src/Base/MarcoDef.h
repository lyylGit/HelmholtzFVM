#pragma once

#define FACEVECTOR(var,V,faceId) Vector<T> var(V[faceId*3],V[faceId*3+1],V[faceId*3+2])