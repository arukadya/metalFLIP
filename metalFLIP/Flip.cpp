//
//  Flip.cpp
//  Mable
//
//  Created by 須之内俊樹 on 2023/02/26.
//

#include "Flip.hpp"

void PIC_FLIP::execute(std::string foldername,std::string filename){
    int cnt = 0;
    for(unsigned int i=0;i<repeatCount;i++){
        if(timer == 2)TD.startTimer("execute");
        if(timer == 1)TD.startTimer("preprocess");
        
        preprocessingParticles();
        if(timer == 1)TD.endTimer();
        particlesVelocityToGrid();
        if(timer == 1)TD.startTimer("GridPressure");
        calGridPressure();
        if(timer == 1)TD.endTimer();
        gridVelocityToParticles();
        advectParticles();
        
        if(i%(repeatCount/500) == 0){
            if(timer == 2){
                TD.endTimer();
                TD.startTimer("output");
            }
            output(vertices);
            implicit_function = cal_implicitFunction(particles, map, radius, dx, Nx, Ny, Nz);
            std::vector<double>data = implicit_function.convert2Vector();
            ImplicitFunction<double> imp = ImplicitFunction<double>(Nx*Ny*Nz, Nx, Ny, Nz, dx, dx, dx, data);
            std::string OutputVTK_imp = foldername + "/output" + "/output"+std::to_string(cnt)+".vtk";
            std::string OutputVTK_iso = foldername + "/isosurface"+ "/output"+std::to_string(cnt)+".off";
            outputVTK(OutputVTK_imp.c_str(),implicit_function,dx,Nx,Ny,Nz);
            marching_cubes(marchingVertices, marchingFaces, origin, dist, imp, threshold);
            outputOFF(OutputVTK_iso.c_str(), marchingVertices, marchingFaces);
            
            volumes.push_back(cal_volume(particles, map, radius, dx, Nx, Ny, Nz, threshold,i));
            
            cnt++;
            if(timer == 2)TD.endTimer();
        }
    }
    outputVolume(filename.c_str(),volumes);
    //std::cout << "dx:" << dx << " dt:" << dt << " Nx*Ny*Nz:" << Nx*Ny*Nz << " repeat:" << repeatCount << std::endl;
    std::cout << " gamma:" << gamma << " radius:" << radius << " threshold:" << threshold <<std::endl;
    //std::cout << "extend:" << extend << std::endl;
    std::cout << std::endl;
}
void PIC_FLIP::output(std::vector<Eigen::Vector3d> &v){
    v.resize(particles.size());
    for(int i=0;i<v.size();i++)v[i] = particles[i].position;
}
void PIC_FLIP::initParticles(){
    for(unsigned int i=0;i<Nx;i++){
        for(unsigned int j=0;j<Ny;j++){
            for(unsigned int k=0;k<Nz;k++){
//                    if((i>Nx/4 && i < Nx/4 * 3) && (j > Ny/4) && ((k > Nz/4 && k < Nz/4*3))){
                if((i>Nx/2) && (k > Nz/2)){
                    Eigen::Vector3d v0 = {0.0,0.0,0.0};
                    std::vector<Eigen::Vector3d>pos(8);//1グリッドにn^2個置くのが流儀らしい
                    pos[0] = Eigen::Vector3d{(i+0.25)*dx,(j+0.25)*dx,(k+0.25)*dx};
                    pos[1] = Eigen::Vector3d{(i+0.75)*dx,(j+0.25)*dx,(k+0.25)*dx};
                    pos[2] = Eigen::Vector3d{(i+0.25)*dx,(j+0.75)*dx,(k+0.25)*dx};
                    pos[3] = Eigen::Vector3d{(i+0.75)*dx,(j+0.75)*dx,(k+0.25)*dx};
                    pos[4] = Eigen::Vector3d{(i+0.25)*dx,(j+0.25)*dx,(k+0.75)*dx};
                    pos[5] = Eigen::Vector3d{(i+0.75)*dx,(j+0.25)*dx,(k+0.75)*dx};
                    pos[6] = Eigen::Vector3d{(i+0.25)*dx,(j+0.75)*dx,(k+0.75)*dx};
                    pos[7] = Eigen::Vector3d{(i+0.75)*dx,(j+0.75)*dx,(k+0.75)*dx};
                    for(unsigned int l=0;l<8;l++){
                        //if(i==Nx-1 && j == Ny/2 - 1)std::cout << pos[k].x() << "," << pos[k].y() << std::endl;
                        Particle tmp = Particle(v0,pos[l]);
                        particles.push_back(tmp);
                        weights.push_back(-1);
                    }
                }
            }
        }
    }
}

void PIC_FLIP::locateParticlesOnGrid(std::unordered_map<std::vector<int>,std::vector<int>,ArrayHasher<3>>&map){
    std::unordered_map<std::vector<int>,std::vector<int>,ArrayHasher<3>>m;
    for(int i=0;i<particles.size();i++){
        std::vector<int>key = calGridOrKey(particles[i].position);
        if(m.find(key) == m.end()){
            std::vector<int>contain = {i};
            m.emplace(key,contain);
        }
        else{
            m.at(key).push_back(i);
        }
        particles[i].setGridIndex(key[0], key[1], key[2]);
    }
    map=m;
    //std::cout << "end locating" << std::endl;
}
Eigen::Vector3d PIC_FLIP::calCellSize(){
    return Eigen::Vector3d{dx,dx,dx};
}
std::vector<int> PIC_FLIP::calGridOrKey(Eigen::Vector3d v){
    std::vector<int> key;
    Eigen::Vector3d cellSize = calCellSize();
    key.push_back( floor(( (v.x() ) / cellSize.x() ) ));
    key.push_back( floor(( (v.y() ) / cellSize.y() ) ));
    key.push_back( floor(( (v.z() ) / cellSize.z() ) ));
    return key;
}
void PIC_FLIP::particlesVelocityToGrid(){
//----------------速度も初期化した方がいい気がする----------------------------------------------------------
//グリッドの切り方はよく考えるべき．スタッカート格子に則って，u,v別々にやるのが合ってる気がする．ブランチはきれ（戒め）
//-----------------------------------------------------------------------------------------------
        //P2Gの方法は諸説あり。
//-----------------------------------------------------------------------------------------------
    //初期化
    //copyVelocity();
    for(unsigned int i=0;i<Nx+1;i++){
        for(unsigned int j=0;j<Ny;j++){
            for(unsigned int k=0;k<Nz;k++){
                umi.value[i][j][k] = 0;
                old_u.value[i][j][k] = 0;
                u.value[i][j][k] = 0;
            }
        }
    }
    for(unsigned int i=0;i<Nx;i++){
        for(unsigned int j=0;j<Ny+1;j++){
            for(unsigned int k=0;k<Nz;k++){
                vmi.value[i][j][k] = 0;
                old_v.value[i][j][k] = 0;
                v.value[i][j][k] = 0;
            }
        }
    }
    for(unsigned int i=0;i<Nx;i++){
        for(unsigned int j=0;j<Ny;j++){
            for(unsigned int k=0;k<Nz+1;k++){
                wmi.value[i][j][k] = 0;
                old_w.value[i][j][k] = 0;
                w.value[i][j][k] = 0;
            }
        }
    }
    //u
    for(int i=1;i<Nx+1;i++){
        for(int j=0;j<Ny;j++){
            for(int k=0;k<Nz;k++){
                std::vector<Eigen::Vector3d>gx_list = {{(i-1)*dx,(j+0.5)*dx,(k+0.5)*dx},{(i)*dx,(j+0.5)*dx,(k+0.5)*dx}};
                std::vector<std::vector<int>>key_list = {{i-1,j,k},{i,j,k}};
                bool flg = false;
                if(extend){
                    if(i != 1){
                        gx_list.push_back({(i-2)*dx,(j+0.5)*dx,(k+0.5)*dx});
                        key_list.push_back({i-2,j,k});
                    }
                    if(i != Nx){
                        gx_list.push_back({(i+1)*dx,(j+0.5)*dx,(k+0.5)*dx});
                        key_list.push_back({i+1,j,k});
                    }
                }
                
                for(int l=0;l<gx_list.size();l++){
                    if(map.find(key_list[l]) != map.end()){
                        flg = true;
                        auto val = map.at(key_list[l]);
                        for(auto x:val){
                            Eigen::Vector3d px = particles[x].position;
                            Eigen::Vector3d pv = particles[x].velocity;
                            double weight = weightFunction(px, gx_list[l], dx);
                            umi.value[i][j][k] += weight*mp;
                            u.value[i][j][k] += weight*mp*pv.x();
                        }
                    }
                }
                if(flg){
                    u.value[i][j][k] /= umi.value[i][j][k];
                    old_u.value[i][j][k] = u.value[i][j][k];
                    u.value[i][j][k] += ufi.value[i][j][k] * dt/umi.value[i][j][k];
                }
            }
        }
    }
    //v
    for(int i=0;i<Nx;i++){
        for(int j=1;j<Ny+1;j++){
            for(int k=0;k<Nz;k++){
                std::vector<Eigen::Vector3d>gx_list = {{(i+0.5)*dx,(j-1)*dx,(k+0.5)*dx},{(i+0.5)*dx,(j)*dx,(k+0.5)*dx}};
                std::vector<std::vector<int>>key_list = {{i,j-1,k},{i,j,k}};
                bool flg = false;
                if(extend){
                    if(j != 1){
                        gx_list.push_back({(i+0.5)*dx,(j-2)*dx,(k+0.5)*dx});
                        key_list.push_back({i,j-2,k});
                    }
                    if(j != Ny){
                        gx_list.push_back({(i+0.5)*dx,(j+1)*dx,(k+0.5)*dx});
                        key_list.push_back({i,j+1,k});
                    }
                }
                for(int l=0;l<gx_list.size();l++){
                    if(map.find(key_list[l]) != map.end()){
                        flg = true;
                        auto val = map.at(key_list[l]);
                        for(auto x:val){
                            Eigen::Vector3d px = particles[x].position;
                            Eigen::Vector3d pv = particles[x].velocity;
                            double weight = weightFunction(px, gx_list[l], dx);
                            vmi.value[i][j][k] += weight*mp;
                            v.value[i][j][k] += weight*mp*pv.y();
                        }
                    }
                }
                if(flg){
                    v.value[i][j][k] /= vmi.value[i][j][k];
                    old_v.value[i][j][k] = v.value[i][j][k];
                    v.value[i][j][k] += vfi.value[i][j][k] * dt/vmi.value[i][j][k];
                }
            }
        }
    }
    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            for(int k=1;k<Nz+1;k++){
                std::vector<Eigen::Vector3d>gx_list = {{(i+0.5)*dx,(j+0.5)*dx,(k-1)*dx},{(i+0.5)*dx,(j+0.5)*dx,(k)*dx}};
                std::vector<std::vector<int>>key_list = {{i,j,k-1},{i,j,k}};
                bool flg = false;
                if(extend){
                    if(k != 1){
                        gx_list.push_back({(i+0.5)*dx,(j+0.5)*dx,(k-2)*dx});
                        key_list.push_back({i,j,k-2});
                    }
                    if(k != Ny){
                        gx_list.push_back({(i+0.5)*dx,(j+0.5)*dx,(k+1)*dx});
                        key_list.push_back({i,j,k+1});
                    }
                }
                for(int l=0;l<gx_list.size();l++){
                    if(map.find(key_list[l]) != map.end()){
                        flg = true;
                        auto val = map.at(key_list[l]);
                        for(auto x:val){
                            Eigen::Vector3d px = particles[x].position;
                            Eigen::Vector3d pv = particles[x].velocity;
                            double weight = weightFunction(px, gx_list[l], dx);
                            wmi.value[i][j][k] += weight*mp;
                            w.value[i][j][k] += weight*mp*pv.z();
                        }
                    }
                }
                if(flg){
                    w.value[i][j][k] /= wmi.value[i][j][k];
                    old_w.value[i][j][k] = w.value[i][j][k];
                    w.value[i][j][k] += wfi.value[i][j][k] * dt/wmi.value[i][j][k];
                }
            }
        }
    }
    //v.print();
}

void PIC_FLIP::calGridPressure(){
    //copyVelocity();
    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            for(int k=0;k<Nz;k++){
                std::vector<int>key = {i,j,k};
                if(map.find(key) == map.end()){
                    p.value[i][j][k] = 0;
                }
            }
        }
    }
    project(map);
}

void PIC_FLIP::gridVelocityToParticles(){
//-----------------------これも速度を初期化した方がいい気がする-----------------------------------------------
    for(auto &p:particles){
        p.PIC_velocity = {0,0,0};
        p.FLIP_velocity = p.velocity;
    }
    for(int i=1;i<Nx+1;i++){
        for(int j=0;j<Ny;j++){
            for(int k=0;k<Nz;k++){
                std::vector<Eigen::Vector3d>gx_list = {{(i-1)*dx,(j+0.5)*dx,(k+0.5)*dx},{(i)*dx,(j+0.5)*dx,(k+0.5)*dx}};
                std::vector<std::vector<int>>key_list = {{i-1,j,k},{i,j,k}};
                if(extend){
                    if(i != 1){
                        gx_list.insert(gx_list.begin(),{(i-2)*dx,(j+0.5)*dx,(k+0.5)*dx});
                        key_list.insert(key_list.begin(),{i-2,j,k});
                    }
                    if(i != Nx){
                        gx_list.push_back({(i+1)*dx,(j+0.5)*dx,(k+0.5)*dx});
                        key_list.push_back({i+1,j,k});
                    }
                }
                for(int l=0;l<gx_list.size();l++){
                    if(map.find(key_list[l]) != map.end()){
                        auto val = map.at(key_list[l]);
                        for(auto x:val){
                            Eigen::Vector3d px = particles[x].position;
                            double weight = weightFunction(px, gx_list[l], dx);
                            if(i == 1){
                                particles[x].PIC_velocity.x()+= weight*mp*u.value[i-1+l][j][k];
                                particles[x].FLIP_velocity.x()+= weight*mp*(u.value[i-1+l][j][k]-old_u.value[i-1+l][j][k]);
                            }
                            else{
                                particles[x].PIC_velocity.x()+= weight*mp*u.value[i-2+l][j][k];
                                particles[x].FLIP_velocity.x()+= weight*mp*(u.value[i-2+l][j][k]-old_u.value[i-2+l][j][k]);
                            }
                        }
                    }
                }
            }
        }
    }
    for(int i=0;i<Nx;i++){
        for(int j=1;j<Ny+1;j++){
            for(int k=0;k<Nz;k++){
                std::vector<Eigen::Vector3d>gx_list = {{(i+0.5)*dx,(j-1)*dx,(k+0.5)*dx},{(i+0.5)*dx,(j)*dx,(k+0.5)*dx}};
                std::vector<std::vector<int>>key_list = {{i,j-1,k},{i,j,k}};
                if(extend){
                    if(j != 1){
                        gx_list.insert(gx_list.begin(),{(i+0.5)*dx,(j-2)*dx,(k+0.5)*dx});
                        key_list.insert(key_list.begin(),{i,j-2,k});
                    }
                    if(j != Ny){
                        gx_list.push_back({(i+0.5)*dx,(j+1)*dx,(k+0.5)*dx});
                        key_list.push_back({i,j+1,k});
                    }
                }
                for(int l=0;l<gx_list.size();l++){
                    if(map.find(key_list[l]) != map.end()){
                        auto val = map.at(key_list[l]);
                        for(auto x:val){
                            Eigen::Vector3d px = particles[x].position;
                            double weight = weightFunction(px, gx_list[l], dx);
                            if(j == 1){
                                particles[x].PIC_velocity.y()+= weight*mp*v.value[i][j-1+l][k];
                                particles[x].FLIP_velocity.y()+= weight*mp*(v.value[i][j-1+l][k]-old_v.value[i][j-1+l][k]);
                            }
                            else{
                                particles[x].PIC_velocity.y()+= weight*mp*v.value[i][j-2+l][k];
                                particles[x].FLIP_velocity.y()+= weight*mp*(v.value[i][j-2+l][k]-old_v.value[i][j-2+l][k]);
                            }
                        }
                    }
                }
            }
        }
    }
    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            for(int k=1;k<Nz+1;k++){
                std::vector<Eigen::Vector3d>gx_list = {{(i+0.5)*dx,(j+0.5)*dx,(k-1)*dx},{(i+0.5)*dx,(j+0.5)*dx,(k)*dx}};
                std::vector<std::vector<int>>key_list = {{i,j,k-1},{i,j,k}};
                if(extend){
                    if(k != 1){
                        gx_list.insert(gx_list.begin(),{(i+0.5)*dx,(j+0.5)*dx,(k-2)*dx});
                        key_list.insert(key_list.begin(),{i,j,k-2});
                    }
                    if(k != Ny){
                        gx_list.push_back({(i+0.5)*dx,(j+0.5)*dx,(k+1)*dx});
                        key_list.push_back({i,j,k+1});
                    }
                }
                for(int l=0;l<gx_list.size();l++){
                    if(map.find(key_list[l]) != map.end()){
                        auto val = map.at(key_list[l]);
                        for(auto x:val){
                            Eigen::Vector3d px = particles[x].position;
                            double weight = weightFunction(px, gx_list[l], dx);
                            if(k == 1){
                                particles[x].PIC_velocity.z()+= weight*mp*w.value[i][j][k-1+l];
                                particles[x].FLIP_velocity.z()+= weight*mp*(w.value[i][j][k-1+l]-old_w.value[i][j][k-1+l]);
                            }
                            else{
                                particles[x].PIC_velocity.z()+= weight*mp*w.value[i][j][k-2+l];
                                particles[x].FLIP_velocity.z()+= weight*mp*(w.value[i][j][k-2+l]-old_w.value[i][j][k-2+l]);
                            }
                        }
                    }
                }
            }
        }
    }
//        for(auto &p:particles){
//            std::cout << p.PIC_velocity << std::endl;
//        }
}
void PIC_FLIP::advectParticles(){
    //std::cout << L << std::endl;
    for(int i=0;i<particles.size();i++){
        //std::cout <<"x0:"<< particles[i].position.x() << " " << particles[i].position.y() << std::endl;
        particles[i].velocity = alpha*particles[i].FLIP_velocity + (1 - alpha)*particles[i].PIC_velocity;
        //if(particles[i].velocity.norm() * dt > dx)std::cout << "CFLError" << std::endl;
        particles[i].position.x() += particles[i].velocity.x()*dt;
        particles[i].position.y() += particles[i].velocity.y()*dt;
        particles[i].position.z() += particles[i].velocity.z()*dt;
        pushout(particles[i].position, L,dx);
        //std::cout << particles[i].velocity.x() << " " << particles[i].velocity.y() << std::endl;
    }
}
void PIC_FLIP::preprocessingParticles(){
    locateParticlesOnGrid(map);
    //各粒子について，４近傍の密度分布を修正するベクトルの計算
    for(int i=0;i<particles.size();i++){
        std::vector<int>key = {particles[i].gridIndex[0],particles[i].gridIndex[1],particles[i].gridIndex[2]};
        std::vector<bool>F = {
            key[0]<Nx-1,
            key[1]<Ny-1,
            key[0]>0,
            key[1]>0,
            key[2]>0,
            key[2]<Nz-1};
        std::vector<std::vector<int>>keys ={
            {key[0]+1,key[1],key[2]},
            {key[0],key[1]+1,key[2]},
            {key[0]-1,key[1],key[2]},
            {key[0],key[1]-1,key[2]},
            {key[0],key[1],key[2]-1},
            {key[0],key[1],key[2]+1}};
        auto val = map.at(key);
        particles[i].fixVector = {0,0,0};
        for(auto &j:val){
            if(j == i)continue;
            particles[i].fixVector += fixDensityVector(particles[j].position, particles[i].position, radius, gamma, dt);
        }
        for(int k=0;k<keys.size();k++){
            if(F[k] && map.find(keys[k]) != map.end()){
                auto val = map.at(keys[k]);
                //particles[i].fixVector = {0,0};
                for(auto &j:val){
                    if(j == i)continue;
                    particles[i].fixVector += fixDensityVector(particles[j].position, particles[i].position, radius, gamma, dt);
                }
            }
        }
//            std::cout << particles[i].fixVector.x() << "," << particles[i].fixVector.y() << std::endl;
    }
    //再サンプリングして速度を修正
    for(int i=0;i<particles.size();i++){
        particles[i].fixVelocity = {0,0,0};
        std::vector<int>key = {particles[i].gridIndex[0],particles[i].gridIndex[1],particles[i].gridIndex[2]};
        std::vector<bool>F = {
            key[0]<Nx-1,
            key[1]<Ny-1,
            key[0]>0,
            key[1]>0,
            key[2]>0,
            key[2]<Nz-1};
        std::vector<std::vector<int>>keys = {
            {key[0]+1,key[1],key[2]},
            {key[0],key[1]+1,key[2]},
            {key[0]-1,key[1],key[2]},
            {key[0],key[1]-1,key[2]},
            {key[0],key[1],key[2]-1},
            {key[0],key[1],key[2]+1}};
        auto val = map.at(key);
        Eigen::Vector3d sumA = {0,0,0};
        double sumB = 0;
        //bool flg = false;
        for(auto &j:val){
            if(j == i)continue;
            //std::cout << particles[i].fixVector.x() << "," << particles[i].fixVector.y() << std::endl;
            //fixParticleVelocity(particles[i], particles[j], radius, sumA, sumB, mp);
            sumB += fixParticleVelocity(particles[i], particles[j], radius, mp);
            sumA += fixParticleVelocity(particles[i], particles[j], radius, mp)*particles[j].velocity;
            //flg = true;
        }
        for(int k=0;k<keys.size();k++){
            if(F[k] && map.find(keys[k]) != map.end()){
                auto val = map.at(keys[k]);
                for(auto &j:val){
                    if(j == i)continue;
                    //fixParticleVelocity(particles[i], particles[j], radius, sumA, sumB, mp);
                    sumB += fixParticleVelocity(particles[i], particles[j], radius, mp);
                    sumA += fixParticleVelocity(particles[i], particles[j], radius, mp)*particles[j].velocity;
                    //flg = true;
                    //std::cout << particles[j].velocity.x() << "," << particles[i].velocity.y()<< std::endl;
                }
            }
        }
        if(sumB > 1.0e-4){
            //std::cout << sumA.x() << "," << sumA.y() << "," << sumB << std::endl;
            particles[i].fixVelocity = sumA/sumB;
            particles[i].velocity = particles[i].fixVelocity;
        }
    }
    //計算したベクトルで位置を修正
    for(int i=0;i<particles.size();i++){
        particles[i].position += particles[i].fixVector;
        pushout(particles[i].position, L,dx);
    }
    locateParticlesOnGrid(map);
}

