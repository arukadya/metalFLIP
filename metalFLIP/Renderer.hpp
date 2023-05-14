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
#include "math.hpp"

static constexpr size_t kInstanceRows = 10;
static constexpr size_t kInstanceColumns = 10;
static constexpr size_t kInstanceDepth = 10;
static constexpr size_t kNumInstances = (kInstanceRows * kInstanceColumns * kInstanceDepth);
static constexpr size_t kMaxFramesInFlight = 3;

class Renderer
{
    public:
        Renderer( MTL::Device* pDevice );
        ~Renderer();
        void buildShaders();
        void buildComputePipeline();
        void buildDepthStencilStates();
        void buildBuffers();
        void draw( MTK::View* pView );
        double dx = 0.1;//セルの大きさ
        double dt = 0.01;//時間の刻み幅
        double rho = 1.0;
        double radius = dx;
        double thr = 1.0;
        double gamma = 1.5;
    
    private:
        float _angle;
        int _frame;
        dispatch_semaphore_t _semaphore;
        dispatch_semaphore_t _sema;
        static const int kMaxFramesInFlight = 3;
        uint _animationIndex;
        MTL::Device* _pDevice;
        MTL::CommandQueue* _pCommandQueue;
        MTL::Library* _pShaderLibrary;
        MTL::RenderPipelineState* _pPSO;
        MTL::ComputePipelineState* _pComputePSO;
        MTL::DepthStencilState* _pDepthStencilState;
        MTL::Buffer* _pVertexDataBuffer;
        MTL::Buffer* _pInstanceDataBuffer[kMaxFramesInFlight];
        MTL::Buffer* _pCameraDataBuffer[kMaxFramesInFlight];
        MTL::Buffer* _pIndexBuffer;
        PIC_FLIP simulator;
};
namespace shader_types
{
    struct VertexData
    {
        simd::float3 position;
        simd::float3 normal;
        simd::float2 texcoord;
    };

    struct InstanceData
    {
        simd::float4x4 instanceTransform;
        simd::float3x3 instanceNormalTransform;
        simd::float4 instanceColor;
    };

    struct CameraData
    {
        simd::float4x4 perspectiveTransform;
        simd::float4x4 worldTransform;
        simd::float3x3 worldNormalTransform;
    };
}

#endif /* Renderer_hpp */
