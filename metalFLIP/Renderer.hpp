//
//  Renderer.hpp
//  metalFLIP
//
//  Created by 須之内俊樹 on 2023/03/23.
//

#ifndef Renderer_hpp
#define Renderer_hpp

#include <Metal/Metal.hpp>
#include <AppKit/AppKit.hpp>
#include <MetalKit/MetalKit.hpp>

#include <simd/simd.h>
#include "Flip.hpp"
#include <sys/stat.h>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include "marching_cubes.hpp"
#include "ioVTK.hpp"
#include "ioOFF.hpp"

static constexpr size_t kNumInstances = 16;
static constexpr size_t kMaxFramesInFlight = 3;

class Renderer
{
    public:
        Renderer( MTL::Device* pDevice );
        ~Renderer();
        void buildShaders();
        void buildBuffers();
        void draw( MTK::View* pView );

    private:
        MTL::Device* _pDevice;
        MTL::CommandQueue* _pCommandQueue;
        MTL::Library* _pShaderLibrary;
        MTL::RenderPipelineState* _pPSO;
        MTL::Buffer* _pVertexDataBuffer;
        MTL::Buffer* _pInstanceDataBuffer[kMaxFramesInFlight];
        MTL::Buffer* _pIndexBuffer;
        float _angle;
        int _frame;
        dispatch_semaphore_t _semaphore;
        static const int kMaxFramesInFlight;
};

#endif /* Renderer_hpp */
