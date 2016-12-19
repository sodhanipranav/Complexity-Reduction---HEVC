/* The copyright in this software is being made available under the BSD
 * License, included below. This software may be subject to other third party
 * and contributor rights, including patent rights, and no such rights are
 * granted under this license.  
 *
 * Copyright (c) 2010-2014, ITU/ISO/IEC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *  * Neither the name of the ITU/ISO/IEC nor the names of its contributors may
 *    be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

/** \file     TEncCu.cpp
    \brief    Coding Unit (CU) encoder class
*/



#include <stdio.h>
#include <time.h>
#include "TEncTop.h"
#include "TEncCu.h"
#include "TEncAnalyze.h"
#include <string.h>


#include <cmath>
#include <algorithm>
using namespace std;
char buffer[50], buffer1[50];
int input1=99;
int c1 = 0, c2 = 0, c3 = 0, c4=0, c5=0;
int d1 = 0, d2=0, d3=0;
int signal1=0, signal2=0;
int check1=0;
double residualy[480][840], residualy_new[480][840], ctucost[100], zerodef[480][840];
int counter=0;
int flag1= 0, flag2= 0, flag3= 0, flag4 =0, numCTU=0;
int horAMP_2NN = 0, verAMP_2NN = 0, noAMP_2NN= 0;
int horAMP_N2N= 0, verAMP_N2N= 0, noAMP_N2N= 0;
int horAMP_NN= 0, verAMP_NN= 0, noAMP_NN= 0;
int horAMP_amp= 0, verAMP_amp= 0, noAMP_amp= 0;
int hor[12]={}, ver[12]={}, no[12]={};
float mtemp;
int fnum, ftum=0, prev=0, pres=3,framecount=0;

//! \ingroup TLibEncoder
//! \{

// ====================================================================================================================
// Constructor / destructor / create / destroy
// ====================================================================================================================

/**
 \param    uiTotalDepth  total number of allowable depth
 \param    uiMaxWidth    largest CU width
 \param    uiMaxHeight   largest CU height
 */
Void TEncCu::create(UChar uhTotalDepth, UInt uiMaxWidth, UInt uiMaxHeight)
{
  Int i;
  
  m_uhTotalDepth   = uhTotalDepth + 1;
  m_ppcBestCU      = new TComDataCU*[m_uhTotalDepth-1];
  m_ppcTempCU      = new TComDataCU*[m_uhTotalDepth-1];
    
  m_ppcPredYuvBest = new TComYuv*[m_uhTotalDepth-1];
  m_ppcResiYuvBest = new TComYuv*[m_uhTotalDepth-1];
  m_ppcRecoYuvBest = new TComYuv*[m_uhTotalDepth-1];
  m_ppcPredYuvTemp = new TComYuv*[m_uhTotalDepth-1];
  m_ppcResiYuvTemp = new TComYuv*[m_uhTotalDepth-1];
  m_ppcRecoYuvTemp = new TComYuv*[m_uhTotalDepth-1];
  m_ppcOrigYuv     = new TComYuv*[m_uhTotalDepth-1];
  
  UInt uiNumPartitions;
  for( i=0 ; i<m_uhTotalDepth-1 ; i++)
  {
    uiNumPartitions = 1<<( ( m_uhTotalDepth - i - 1 )<<1 );
    UInt uiWidth  = uiMaxWidth  >> i;
    UInt uiHeight = uiMaxHeight >> i;
    
    m_ppcBestCU[i] = new TComDataCU; m_ppcBestCU[i]->create( uiNumPartitions, uiWidth, uiHeight, false, uiMaxWidth >> (m_uhTotalDepth - 1) );
    m_ppcTempCU[i] = new TComDataCU; m_ppcTempCU[i]->create( uiNumPartitions, uiWidth, uiHeight, false, uiMaxWidth >> (m_uhTotalDepth - 1) );
    
    m_ppcPredYuvBest[i] = new TComYuv; m_ppcPredYuvBest[i]->create(uiWidth, uiHeight);
    m_ppcResiYuvBest[i] = new TComYuv; m_ppcResiYuvBest[i]->create(uiWidth, uiHeight);
    m_ppcRecoYuvBest[i] = new TComYuv; m_ppcRecoYuvBest[i]->create(uiWidth, uiHeight);
    
    m_ppcPredYuvTemp[i] = new TComYuv; m_ppcPredYuvTemp[i]->create(uiWidth, uiHeight);
    m_ppcResiYuvTemp[i] = new TComYuv; m_ppcResiYuvTemp[i]->create(uiWidth, uiHeight);
    m_ppcRecoYuvTemp[i] = new TComYuv; m_ppcRecoYuvTemp[i]->create(uiWidth, uiHeight);
    
    m_ppcOrigYuv    [i] = new TComYuv; m_ppcOrigYuv    [i]->create(uiWidth, uiHeight);
  }
  
  m_bEncodeDQP = false;

  // initialize partition order.
  UInt* piTmp = &g_auiZscanToRaster[0];
  initZscanToRaster( m_uhTotalDepth, 1, 0, piTmp);
  initRasterToZscan( uiMaxWidth, uiMaxHeight, m_uhTotalDepth );
  
  // initialize conversion matrix from partition index to pel
  initRasterToPelXY( uiMaxWidth, uiMaxHeight, m_uhTotalDepth );
}

Void TEncCu::destroy()
{
  Int i;
  
  for( i=0 ; i<m_uhTotalDepth-1 ; i++)
  {
    if(m_ppcBestCU[i])
    {
      m_ppcBestCU[i]->destroy();      delete m_ppcBestCU[i];      m_ppcBestCU[i] = NULL;
    }
    if(m_ppcTempCU[i])
    {
      m_ppcTempCU[i]->destroy();      delete m_ppcTempCU[i];      m_ppcTempCU[i] = NULL;
    }
    if(m_ppcPredYuvBest[i])
    {
      m_ppcPredYuvBest[i]->destroy(); delete m_ppcPredYuvBest[i]; m_ppcPredYuvBest[i] = NULL;
    }
    if(m_ppcResiYuvBest[i])
    {
      m_ppcResiYuvBest[i]->destroy(); delete m_ppcResiYuvBest[i]; m_ppcResiYuvBest[i] = NULL;
    }
    if(m_ppcRecoYuvBest[i])
    {
      m_ppcRecoYuvBest[i]->destroy(); delete m_ppcRecoYuvBest[i]; m_ppcRecoYuvBest[i] = NULL;
    }
    if(m_ppcPredYuvTemp[i])
    {
      m_ppcPredYuvTemp[i]->destroy(); delete m_ppcPredYuvTemp[i]; m_ppcPredYuvTemp[i] = NULL;
    }
    if(m_ppcResiYuvTemp[i])
    {
      m_ppcResiYuvTemp[i]->destroy(); delete m_ppcResiYuvTemp[i]; m_ppcResiYuvTemp[i] = NULL;
    }
    if(m_ppcRecoYuvTemp[i])
    {
      m_ppcRecoYuvTemp[i]->destroy(); delete m_ppcRecoYuvTemp[i]; m_ppcRecoYuvTemp[i] = NULL;
    }
    if(m_ppcOrigYuv[i])
    {
      m_ppcOrigYuv[i]->destroy();     delete m_ppcOrigYuv[i];     m_ppcOrigYuv[i] = NULL;
    }
  }
  if(m_ppcBestCU)
  {
    delete [] m_ppcBestCU;
    m_ppcBestCU = NULL;
  }
  if(m_ppcTempCU)
  {
    delete [] m_ppcTempCU;
    m_ppcTempCU = NULL;
  }
  
  if(m_ppcPredYuvBest)
  {
    delete [] m_ppcPredYuvBest;
    m_ppcPredYuvBest = NULL;
  }
  if(m_ppcResiYuvBest)
  {
    delete [] m_ppcResiYuvBest;
    m_ppcResiYuvBest = NULL;
  }
  if(m_ppcRecoYuvBest)
  {
    delete [] m_ppcRecoYuvBest;
    m_ppcRecoYuvBest = NULL;
  }
  if(m_ppcPredYuvTemp)
  {
    delete [] m_ppcPredYuvTemp;
    m_ppcPredYuvTemp = NULL;
  }
  if(m_ppcResiYuvTemp)
  {
    delete [] m_ppcResiYuvTemp;
    m_ppcResiYuvTemp = NULL;
  }
  if(m_ppcRecoYuvTemp)
  {
    delete [] m_ppcRecoYuvTemp;
    m_ppcRecoYuvTemp = NULL;
  }
  if(m_ppcOrigYuv)
  {
    delete [] m_ppcOrigYuv;
    m_ppcOrigYuv = NULL;
  }
}

/** \param    pcEncTop      pointer of encoder class
 */
Void TEncCu::init( TEncTop* pcEncTop )
{
  m_pcEncCfg           = pcEncTop;
  m_pcPredSearch       = pcEncTop->getPredSearch();
  m_pcTrQuant          = pcEncTop->getTrQuant();
  m_pcBitCounter       = pcEncTop->getBitCounter();
  m_pcRdCost           = pcEncTop->getRdCost();
  
  m_pcEntropyCoder     = pcEncTop->getEntropyCoder();
  m_pcCavlcCoder       = pcEncTop->getCavlcCoder();
  m_pcSbacCoder       = pcEncTop->getSbacCoder();
  m_pcBinCABAC         = pcEncTop->getBinCABAC();
  
  m_pppcRDSbacCoder   = pcEncTop->getRDSbacCoder();
  m_pcRDGoOnSbacCoder = pcEncTop->getRDGoOnSbacCoder();
  
  m_pcRateCtrl        = pcEncTop->getRateCtrl();
}

// ====================================================================================================================
// Public member functions
// ====================================================================================================================

/** \param  rpcCU pointer of CU data class
 */
Void TEncCu::compressCU( TComDataCU*& rpcCU )
{
  // initialize CU data
  m_ppcBestCU[0]->initCU( rpcCU->getPic(), rpcCU->getAddr() );
  m_ppcTempCU[0]->initCU( rpcCU->getPic(), rpcCU->getAddr() );

  // analysis of CU
  xCompressCU( m_ppcBestCU[0], m_ppcTempCU[0], 0 );

#if ADAPTIVE_QP_SELECTION
  if( m_pcEncCfg->getUseAdaptQpSelect() )
  {
    if(rpcCU->getSlice()->getSliceType()!=I_SLICE) //IIII
    {
      xLcuCollectARLStats( rpcCU);
    }
  }
#endif
}
/** \param  pcCU  pointer of CU data class
 */
Void TEncCu::encodeCU ( TComDataCU* pcCU )
{
  if ( pcCU->getSlice()->getPPS()->getUseDQP() )
  {
    setdQPFlag(true);
  }

  // Encode CU data
  xEncodeCU( pcCU, 0, 0 );
}

// ====================================================================================================================
// Protected member functions
// ====================================================================================================================
/** Derive small set of test modes for AMP encoder speed-up
 *\param   rpcBestCU
 *\param   eParentPartSize
 *\param   bTestAMP_Hor
 *\param   bTestAMP_Ver
 *\param   bTestMergeAMP_Hor
 *\param   bTestMergeAMP_Ver
 *\returns Void 
*/


#if AMP_ENC_SPEEDUP
#if AMP_MRG
Void TEncCu::deriveTestModeAMP (TComDataCU *&rpcBestCU, PartSize eParentPartSize, Bool &bTestAMP_Hor, Bool &bTestAMP_Ver, Bool &bTestMergeAMP_Hor, Bool &bTestMergeAMP_Ver)
#else
Void TEncCu::deriveTestModeAMP (TComDataCU *&rpcBestCU, PartSize eParentPartSize, Bool &bTestAMP_Hor, Bool &bTestAMP_Ver)
#endif
{
    
    if ((eParentPartSize==SIZE_2NxN) && (rpcBestCU->getPartitionSize(0)==SIZE_2NxN))
    {
        bTestAMP_Hor = true;
    }
    
    else if ((eParentPartSize==SIZE_2NxN) && (rpcBestCU->getPartitionSize(0) == SIZE_Nx2N))
    {
        
        bTestAMP_Ver = true;
        
    }
    
    else if ((eParentPartSize==SIZE_2NxN) && (rpcBestCU->getPartitionSize(0) == SIZE_2Nx2N))
    {
        
        
    }
    
    else if ((eParentPartSize==SIZE_Nx2N) && (rpcBestCU->getPartitionSize(0) == SIZE_2NxN))
    {
        bTestAMP_Hor = true;
        //bTestMergeAMP_Ver = true;
        
    }
    
    else if ((eParentPartSize==SIZE_Nx2N) && (rpcBestCU->getPartitionSize(0) == SIZE_Nx2N))
    {
        bTestAMP_Ver = true;
        
    }
    
    else if ((eParentPartSize==SIZE_Nx2N) && (rpcBestCU->getPartitionSize(0) == SIZE_2Nx2N))
    {
        
        //bTestMergeAMP_Ver = true;
 
        
    }
    
    else if ((eParentPartSize==SIZE_2Nx2N) && (rpcBestCU->getPartitionSize(0) == SIZE_2NxN))
    {
        //bTestAMP_Hor = true;
        //bTestMergeAMP_Ver = true;
    }
    
    else if ((eParentPartSize==SIZE_2Nx2N) && (rpcBestCU->getPartitionSize(0) == SIZE_Nx2N))
    {
        //bTestAMP_Ver = true;
        //bTestMergeAMP_Hor = true;
        
    }
    
    else if ((eParentPartSize==SIZE_2Nx2N) && (rpcBestCU->getPartitionSize(0) == SIZE_2Nx2N))
    {
        
        
    }
    
    else if ((eParentPartSize == SIZE_2NxnU || eParentPartSize == SIZE_nRx2N || eParentPartSize == SIZE_2NxnD || eParentPartSize == SIZE_nLx2N) && (rpcBestCU->getPartitionSize(0) == SIZE_2NxN))
    {
        bTestAMP_Hor = true;
        //bTestMergeAMP_Ver = true;
    }
    
    else if ((eParentPartSize == SIZE_2NxnU || eParentPartSize == SIZE_nRx2N || eParentPartSize == SIZE_2NxnD || eParentPartSize == SIZE_nLx2N) && (rpcBestCU->getPartitionSize(0) == SIZE_Nx2N))
    {
        bTestAMP_Ver = true;
        //bTestMergeAMP_Hor = true;
        
    }
    
    else if ((eParentPartSize == SIZE_2NxnU || eParentPartSize == SIZE_nRx2N || eParentPartSize == SIZE_2NxnD || eParentPartSize == SIZE_nLx2N) && (rpcBestCU->getPartitionSize(0) == SIZE_2Nx2N))
    {
        //bTestMergeAMP_Hor = true;
        //bTestMergeAMP_Ver = true;
    }
    
    if ( rpcBestCU->getWidth(0) == 64 )
    {
        bTestAMP_Hor = false;
        bTestAMP_Ver = false;
    }
  
    /*
  if ( rpcBestCU->getPartitionSize(0) == SIZE_2NxN )
  {
    bTestAMP_Hor = true;
  
  }
  
  else if ( rpcBestCU->getPartitionSize(0) == SIZE_Nx2N )
  {
    bTestAMP_Ver = true;
  
  }
  
  else if ( rpcBestCU->getPartitionSize(0) == SIZE_2Nx2N && rpcBestCU->getMergeFlag(0) == false && rpcBestCU->isSkipped(0) == false )
  {
    bTestAMP_Hor = true;
    bTestAMP_Ver = true;
  
   
  }
    printf("c1 = %d \n", c1);
    printf("c2 = %d \n", c2);
    printf("c3 = %d \n", c3);
    printf("c4 = %d \n", c4);
    printf("c5 = %d \n", c5);
    printf("c6 = %d \n \n", c6);
     */
    
    
#if AMP_MRG
  //! Utilizing the partition size of parent PU
  if ( eParentPartSize >= SIZE_2NxnU && eParentPartSize <= SIZE_nRx2N )
  {
    bTestMergeAMP_Hor = true;
    bTestMergeAMP_Ver = true;
  }

  if ( eParentPartSize == SIZE_NONE ) //! if parent is intra
  {
    if ( rpcBestCU->getPartitionSize(0) == SIZE_2NxN )
    {
      bTestMergeAMP_Hor = true;
    }
    else if ( rpcBestCU->getPartitionSize(0) == SIZE_Nx2N )
    {
      bTestMergeAMP_Ver = true;
    }
  }

  if ( rpcBestCU->getPartitionSize(0) == SIZE_2Nx2N && rpcBestCU->isSkipped(0) == false )
  {
    bTestMergeAMP_Hor = true;
    bTestMergeAMP_Ver = true;
  }

  if ( rpcBestCU->getWidth(0) == 64 )
  {
      //c1++;
      //printf("c1 = %d \n", c1);
      bTestAMP_Hor = false;
      bTestAMP_Ver = false;
  }
    
#else
  /*
   
   //! Utilizing the partition size of parent PU
  if ( eParentPartSize >= SIZE_2NxnU && eParentPartSize <= SIZE_nRx2N )
  {
    bTestAMP_Hor = true;
    bTestAMP_Ver = true;
  }

  if ( eParentPartSize == SIZE_2Nx2N )
  {
    bTestAMP_Hor = false;
    bTestAMP_Ver = false;
  }
   
   */
    
#endif
}
#endif

// ====================================================================================================================
// Protected member functions
// ====================================================================================================================
/** Compress a CU block recursively with enabling sub-LCU-level delta QP
 *\param   rpcBestCU
 *\param   rpcTempCU
 *\param   uiDepth
 *\returns Void
 *
 *- for loop of QP value to compress the current CU with all possible QP
*/
#if AMP_ENC_SPEEDUP
Void TEncCu::xCompressCU( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, UInt uiDepth, PartSize eParentPartSize )
#else
Void TEncCu::xCompressCU( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, UInt uiDepth, PartSize eParentPartSize  )
#endif
{
    TComPic* pcPic = rpcBestCU->getPic();

  // get Original YUV data from picture
  m_ppcOrigYuv[uiDepth]->copyFromPicYuv( pcPic->getPicYuvOrg(), rpcBestCU->getAddr(), rpcBestCU->getZorderIdxInCU() );

  // variable for Early CU determination
  Bool    bSubBranch = true;

  // variable for Cbf fast mode PU decision
  Bool    doNotBlockPu = true;
  Bool earlyDetectionSkipMode = false;

  Bool bBoundary = false;
  UInt uiLPelX   = rpcBestCU->getCUPelX();
  UInt uiRPelX   = uiLPelX + rpcBestCU->getWidth(0)  - 1;
  UInt uiTPelY   = rpcBestCU->getCUPelY();
  UInt uiBPelY   = uiTPelY + rpcBestCU->getHeight(0) - 1;

  Int iBaseQP = xComputeQP( rpcBestCU, uiDepth );
  Int iMinQP;
  Int iMaxQP;
  Bool isAddLowestQP = false;

  if( (g_uiMaxCUWidth>>uiDepth) >= rpcTempCU->getSlice()->getPPS()->getMinCuDQPSize() )
  {
    Int idQP = m_pcEncCfg->getMaxDeltaQP();
    iMinQP = Clip3( -rpcTempCU->getSlice()->getSPS()->getQpBDOffsetY(), MAX_QP, iBaseQP-idQP );
    iMaxQP = Clip3( -rpcTempCU->getSlice()->getSPS()->getQpBDOffsetY(), MAX_QP, iBaseQP+idQP );
  }
  else
  {
    iMinQP = rpcTempCU->getQP(0);
    iMaxQP = rpcTempCU->getQP(0);
  }

  if ( m_pcEncCfg->getUseRateCtrl() )
  {
    iMinQP = m_pcRateCtrl->getRCQP();
    iMaxQP = m_pcRateCtrl->getRCQP();
  }

  // transquant-bypass (TQB) processing loop variable initialisation ---

  const Int lowestQP = iMinQP; // For TQB, use this QP which is the lowest non TQB QP tested (rather than QP'=0) - that way delta QPs are smaller, and TQB can be tested at all CU levels.

  if ( (rpcTempCU->getSlice()->getPPS()->getTransquantBypassEnableFlag()) )
  {
    isAddLowestQP = true; // mark that the first iteration is to cost TQB mode.
    iMinQP = iMinQP - 1;  // increase loop variable range by 1, to allow testing of TQB mode along with other QPs
    if ( m_pcEncCfg->getCUTransquantBypassFlagForceValue() )
    {
      iMaxQP = iMinQP;
    }
  }

  // If slice start or slice end is within this cu...
    
  TComSlice * pcSlice = rpcTempCU->getPic()->getSlice(rpcTempCU->getPic()->getCurrSliceIdx());
    
  Bool bSliceStart = pcSlice->getSliceSegmentCurStartCUAddr()>rpcTempCU->getSCUAddr()&&pcSlice->getSliceSegmentCurStartCUAddr()<rpcTempCU->getSCUAddr()+rpcTempCU->getTotalNumPart();
    
  Bool bSliceEnd = (pcSlice->getSliceSegmentCurEndCUAddr()>rpcTempCU->getSCUAddr()&&pcSlice->getSliceSegmentCurEndCUAddr()<rpcTempCU->getSCUAddr()+rpcTempCU->getTotalNumPart());
    
  Bool bInsidePicture = ( uiRPelX < rpcBestCU->getSlice()->getSPS()->getPicWidthInLumaSamples() ) && ( uiBPelY < rpcBestCU->getSlice()->getSPS()->getPicHeightInLumaSamples() );
    
  // We need to split, so don't try these modes.
    
    
    
  if(!bSliceEnd && !bSliceStart && bInsidePicture )
  {
    for (Int iQP=iMinQP; iQP<=iMaxQP; iQP++)
    {
      const Bool bIsLosslessMode = isAddLowestQP && (iQP == iMinQP);

      if (bIsLosslessMode)
      {
        iQP = lowestQP;
      }

      rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
        
        
      // do inter modes, SKIP and 2Nx2N
        
        
      if( rpcBestCU->getSlice()->getSliceType() != I_SLICE )
      {
        // 2Nx2N
        if(m_pcEncCfg->getUseEarlySkipDetection())
        {
          xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_2Nx2N, uiDepth, signal1 );
          rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );//by Competition for inter_2Nx2N
        }
        // SKIP
        xCheckRDCostMerge2Nx2N( rpcBestCU, rpcTempCU, &earlyDetectionSkipMode );//by Merge for inter_2Nx2N
        rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );

        if(!m_pcEncCfg->getUseEarlySkipDetection())
        {
          // 2Nx2N, NxN
            signal1=0;
            
          xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_2Nx2N, uiDepth, signal1 );
          rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
          if(m_pcEncCfg->getUseCbfFastMode())
          {
            doNotBlockPu = rpcBestCU->getQtRootCbf( 0 ) != 0;
              
          }
        }
    }
        
       if (bIsLosslessMode)
      {
        iQP = iMinQP;
      }
    }
      
      
    if(!earlyDetectionSkipMode)
    {
      for (Int iQP=iMinQP; iQP<=iMaxQP; iQP++)
      {
        const Bool bIsLosslessMode = isAddLowestQP && (iQP == iMinQP);

        if (bIsLosslessMode)
        {
          iQP = lowestQP;
        }
        rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );

        // do inter modes, NxN, 2NxN, and Nx2N
        if( rpcBestCU->getSlice()->getSliceType() != I_SLICE )
        {
          // 2Nx2N, NxN
            
          if(!( (rpcBestCU->getWidth(0)==8) && (rpcBestCU->getHeight(0)==8) ))
          {
            if( uiDepth == g_uiMaxCUDepth - g_uiAddCUDepth && doNotBlockPu)
            {
             xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_NxN, uiDepth, signal1);
             rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
            }
          }

          // 2NxN, Nx2N
          if(doNotBlockPu)
          {
             
            
            xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_Nx2N, uiDepth, signal1);
            rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
            if(m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_Nx2N )
            {
              doNotBlockPu = rpcBestCU->getQtRootCbf( 0 ) != 0;
            }
          }
          if(doNotBlockPu)
          {
              
            xCheckRDCostInter      ( rpcBestCU, rpcTempCU, SIZE_2NxN, uiDepth, signal1);
            rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
            if(m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_2NxN)
            {
              doNotBlockPu = rpcBestCU->getQtRootCbf( 0 ) != 0;
            }
          }
            if(rpcBestCU->getSlice()->getPOC()==pres+1){
                ftum=0;
                pres++;}
            
            if(rpcBestCU->getSlice()->getPOC()==pres && ftum==0)
            {
                c1=0;
                d1=0;
                flag1=0;
                flag2=0;
                signal1=0;
                ::check1=0;
                ::numCTU=0;
                ::counter = 0;
                ftum++;
                
            }

            
            if(!m_pcEncCfg->getUseEarlySkipDetection() && numCTU<=91)
            {
                // 2Nx2N, NxN
                signal1=1;
                xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_2Nx2N, uiDepth, signal1 );
                signal1=0;
                //flag1++;
                //printf("Counter 1 : Checked in %d \n\n", flag1);
                rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
                if(m_pcEncCfg->getUseCbfFastMode())
                {
                    doNotBlockPu = rpcBestCU->getQtRootCbf( 0 ) != 0;
                }
            }
          
            
            #if 1
          //! Try AMP (SIZE_2NxnU, SIZE_2NxnD, SIZE_nLx2N, SIZE_nRx2N)
          if( pcPic->getSlice(0)->getSPS()->getAMPAcc(uiDepth) )
          {
              
#if AMP_ENC_SPEEDUP        
            Bool bTestAMP_Hor = false, bTestAMP_Ver = false;
              

#if AMP_MRG
            Bool bTestMergeAMP_Hor = false, bTestMergeAMP_Ver = false;

            deriveTestModeAMP (rpcBestCU, eParentPartSize, bTestAMP_Hor, bTestAMP_Ver, bTestMergeAMP_Hor, bTestMergeAMP_Ver);
#else
            deriveTestModeAMP (rpcBestCU, eParentPartSize, bTestAMP_Hor, bTestAMP_Ver);
#endif
           
             /* 
              
              if(bTestAMP_Hor && !bTestAMP_Ver)
                  c1++;
              if(bTestAMP_Ver && !bTestAMP_Hor )
                  c2++;
              if(!bTestAMP_Hor && !bTestAMP_Ver)
                  c3++;
              if(bTestAMP_Hor && bTestAMP_Ver)
                  c4++;
              
              mtemp = c1;
              mtemp = 100*mtemp/(c1+c2+c3+c4);
              printf("Hor only = %f \n", mtemp);
              
              mtemp = c2;
              mtemp = 100*mtemp/(c1+c2+c3+c4);
              printf("Ver only = %f \n", mtemp);
              
              mtemp = c3;
              mtemp = 100*mtemp/(c1+c2+c3+c4);
              printf("None = %f \n", mtemp);
              
              mtemp = c4;
              mtemp = 100*mtemp/(c1+c2+c3+c4);
              printf("Both = %f \n \n", mtemp);
              
              //
              
              if(!bTestAMP_Hor && bTestMergeAMP_Hor && !bTestMergeAMP_Ver) // hor merge only
                  c5++;
              if(!bTestAMP_Ver && bTestMergeAMP_Ver && !bTestMergeAMP_Hor ) // ver merge only
                  c6++;
              if((!bTestMergeAMP_Hor && !bTestMergeAMP_Ver) || (bTestAMP_Hor && bTestAMP_Ver)) // neither hor merge nor ver merge
                  c7++;
              if(bTestMergeAMP_Hor && bTestMergeAMP_Ver && !bTestAMP_Hor && !bTestAMP_Ver) // both hor and ver merge
                  c8++;
              
              mtemp = c5;
              mtemp = 100*mtemp/(c5+c6+c7+c8);
              printf("Hor only = %f \n", mtemp);
              
              mtemp = c6;
              mtemp = 100*mtemp/(c5+c6+c7+c8);
              printf("Ver only = %f \n", mtemp);
              
              mtemp = c7;
              mtemp = 100*mtemp/(c5+c6+c7+c8);
              printf("None = %f \n", mtemp);
              
              mtemp = c8;
              mtemp = 100*mtemp/(c5+c6+c7+c8);
              printf("Both = %f \n \n", mtemp);*/
    
            //! Do horizontal AMP
            if ( bTestAMP_Hor )
            {
              if(doNotBlockPu)
              {
                  //clock_t t;
                  //t = clock();
                  xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_2NxnU, signal1 );
                  //t = clock() - t;
                  //c1 = c1 + t;
                  //printf("time taken, Nx2N = %f \n", ((float)c1)/CLOCKS_PER_SEC);
                
                rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
                if(m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_2NxnU )
                {
                  doNotBlockPu = rpcBestCU->getQtRootCbf( 0 ) != 0;
                }
              }
              if(doNotBlockPu)
              {
                  //clock_t t;
                  //t = clock();
                  xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_2NxnD, signal1 );
                  ////t = clock() - t;
                  //c2 = c2 + t;
                  //printf("time taken, Nx2N = %f \n", ((float)c2)/CLOCKS_PER_SEC);
               
                rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
                if(m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_2NxnD )
                {
                  doNotBlockPu = rpcBestCU->getQtRootCbf( 0 ) != 0;
                }
              }
            }
              
              
 #if AMP_MRG
            else if ( bTestMergeAMP_Hor ) 
            {
              if(doNotBlockPu)
              {
                xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_2NxnU, true );
                rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
                if(m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_2NxnU )
                {
                  doNotBlockPu = rpcBestCU->getQtRootCbf( 0 ) != 0;
                }
              }
              if(doNotBlockPu)
              {
                xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_2NxnD, true );
                rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
                if(m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_2NxnD )
                {
                  doNotBlockPu = rpcBestCU->getQtRootCbf( 0 ) != 0;
                }
              }
            }
#endif

            //! Do Vertical AMP
            if ( bTestAMP_Ver )
            {
              if(doNotBlockPu)
              {
                  
                xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_nLx2N, signal1 );
                rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
                if(m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_nLx2N )
                {
                  doNotBlockPu = rpcBestCU->getQtRootCbf( 0 ) != 0;
                }
              }
              if(doNotBlockPu)
              {
               xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_nRx2N, signal1 );
                rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
              }
            }
#if AMP_MRG
            else if ( bTestMergeAMP_Ver )
            {
              if(doNotBlockPu)
              {
               xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_nLx2N, true );
                rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
                if(m_pcEncCfg->getUseCbfFastMode() && rpcBestCU->getPartitionSize(0) == SIZE_nLx2N )
                {
                  doNotBlockPu = rpcBestCU->getQtRootCbf( 0 ) != 0;
                }
              }
              if(doNotBlockPu)
              {
              xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_nRx2N, true );
                rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
              }
            }
#endif

#else
            /*
             
             if ((eParentPartSize==SIZE_2NxN) && (rpcBestCU->getPartitionSize(0)==SIZE_2NxN))
              {
                  counter[0]++;
                  for(int iter=0; iter<12; iter++)
                  {
                      flag[iter]=0;
                  }
                  flag[0]=1;
                  
              }
              
              else if ((eParentPartSize==SIZE_2NxN) && (rpcBestCU->getPartitionSize(0) == SIZE_Nx2N))
              {
                  counter[1]++;
                  for(int iter=0; iter<12; iter++)
                  {
                      flag[iter]=0;
                  }
                  flag[1]=1;
                  
              }
              
              else if ((eParentPartSize==SIZE_2NxN) && (rpcBestCU->getPartitionSize(0) == SIZE_2Nx2N))
              {
                  counter[2]++;
                  for(int iter=0; iter<12; iter++)
                  {
                      flag[iter]=0;
                  }
                  flag[2]=1;
                  
              }
              
              else if ((eParentPartSize==SIZE_Nx2N) && (rpcBestCU->getPartitionSize(0) == SIZE_2NxN))
              {
                  counter[3]++;
                  for(int iter=0; iter<12; iter++)
                  {
                      flag[iter]=0;
                  }
                  flag[3]=1;
                  
              }
              
              else if ((eParentPartSize==SIZE_Nx2N) && (rpcBestCU->getPartitionSize(0) == SIZE_Nx2N))
              {
                  counter[4]++;
                  for(int iter=0; iter<12; iter++)
                  {
                      flag[iter]=0;
                  }
                  flag[4]=1;
                  
              }
              
              else if ((eParentPartSize==SIZE_Nx2N) && (rpcBestCU->getPartitionSize(0) == SIZE_2Nx2N))
              {
                  counter[5]++;
                  for(int iter=0; iter<12; iter++)
                  {
                      flag[iter]=0;
                  }
                  flag[5]=1;
                  
              }
              
              else if ((eParentPartSize==SIZE_2Nx2N) && (rpcBestCU->getPartitionSize(0) == SIZE_2NxN))
              {
                  counter[6]++;
                  for(int iter=0; iter<12; iter++)
                  {
                      flag[iter]=0;
                  }
                  flag[6]=1;
                  
              }
              
              else if ((eParentPartSize==SIZE_2Nx2N) && (rpcBestCU->getPartitionSize(0) == SIZE_Nx2N))
              {
                  counter[7]++;
                  for(int iter=0; iter<12; iter++)
                  {
                      flag[iter]=0;
                  }
                  flag[7]=1;
                  
              }
              
              else if ((eParentPartSize==SIZE_2Nx2N) && (rpcBestCU->getPartitionSize(0) == SIZE_2Nx2N))
              {
                  counter[8]++;
                  for(int iter=0; iter<12; iter++)
                  {
                      flag[iter]=0;
                  }
                  flag[8]=1;
                  
              }
              
              else if ((eParentPartSize == SIZE_2NxnU || eParentPartSize == SIZE_nRx2N || eParentPartSize == SIZE_2NxnD || eParentPartSize == SIZE_nLx2N) && (rpcBestCU->getPartitionSize(0) == SIZE_2NxN))
              {
                  counter[9]++;
                  for(int iter=0; iter<12; iter++)
                  {
                      flag[iter]=0;
                  }
                  flag[9]=1;
                  
              }
              
              else if ((eParentPartSize == SIZE_2NxnU || eParentPartSize == SIZE_nRx2N || eParentPartSize == SIZE_2NxnD || eParentPartSize == SIZE_nLx2N) && (rpcBestCU->getPartitionSize(0) == SIZE_Nx2N))
              {
                  counter[10]++;
                  for(int iter=0; iter<12; iter++)
                  {
                      flag[iter]=0;
                  }
                  flag[10]=1;
                  
              }
              
              else if ((eParentPartSize == SIZE_2NxnU || eParentPartSize == SIZE_nRx2N || eParentPartSize == SIZE_2NxnD || eParentPartSize == SIZE_nLx2N) && (rpcBestCU->getPartitionSize(0) == SIZE_2Nx2N))
              {
                  counter[11]++;
                  for(int iter=0; iter<12; iter++)
                  {
                      flag[iter]=0;
                  }
                  flag[11]=1;
                  
              }
              
              else
              {
                  c1++;
                  for(int iter=0; iter<12; iter++)
                  {
                      flag[iter]=0;
                  }
              }
             
             */

     
            
            //xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_2NxnU, uiDepth );
            //rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode);
            //xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_2NxnD, uiDepth );
            //rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
            //xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_nLx2N, uiDepth );
            //rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
            //xCheckRDCostInter( rpcBestCU, rpcTempCU, SIZE_nRx2N, uiDepth );
            //rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
       
              
             /* 
              
              for(int iter1=0; iter1<12; iter1++)
              {
                  if(flag[iter1]==1)
                  {
                      if ( rpcBestCU->getPartitionSize(0) == SIZE_2NxnU || rpcBestCU->getPartitionSize(0) == SIZE_2NxnD) // implies horizontal AMP
                          
                          hor[iter1]++;
                      
                      else if ( rpcBestCU->getPartitionSize(0) == SIZE_nLx2N || rpcBestCU->getPartitionSize(0) == SIZE_nLx2N) // implies vertical AMP
                          
                          ver[iter1]++;
                      
                      else // implies NOT AMP
                          
                          no[iter1]++;

                  }
              }
              
              
              for(int iter1=0; iter1<12;iter1++)
              {
                  if(counter[iter1]!=0)
                  {
                      mytemp=hor[iter1];//counter[iter1];
                      mytemp=mytemp/counter[iter1];
                  printf("%d. hor- %f \n",iter1+1, mytemp);
                  
                      mytemp=ver[iter1];//counter[iter1];
                      mytemp=mytemp/counter[iter1];
                  printf("%d. ver- %f \n",iter1+1, mytemp);
                  
                      mytemp=no[iter1];//counter[iter1];
                      mytemp=mytemp/counter[iter1];
                  printf("%d. no- %f \n \n",iter1+1, mytemp);
                  }
              }
              
              printf("c1 = %d \n", c1);
              
              */
#endif
          }
#endif
        }
          
          
        // do normal intra modes
        // speedup for inter frames
        if( rpcBestCU->getSlice()->getSliceType() == I_SLICE || 
          rpcBestCU->getCbf( 0, TEXT_LUMA     ) != 0   ||
          rpcBestCU->getCbf( 0, TEXT_CHROMA_U ) != 0   ||
          rpcBestCU->getCbf( 0, TEXT_CHROMA_V ) != 0     ) // avoid very complex intra if it is unlikely
        {
          xCheckRDCostIntra( rpcBestCU, rpcTempCU, SIZE_2Nx2N );
          rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
          if( uiDepth == g_uiMaxCUDepth - g_uiAddCUDepth )
          {
            if( rpcTempCU->getWidth(0) > ( 1 << rpcTempCU->getSlice()->getSPS()->getQuadtreeTULog2MinSize() ) )
            {
              xCheckRDCostIntra( rpcBestCU, rpcTempCU, SIZE_NxN   );
              rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
            }
          }
        }

        // test PCM
        if(pcPic->getSlice(0)->getSPS()->getUsePCM()
          && rpcTempCU->getWidth(0) <= (1<<pcPic->getSlice(0)->getSPS()->getPCMLog2MaxSize())
          && rpcTempCU->getWidth(0) >= (1<<pcPic->getSlice(0)->getSPS()->getPCMLog2MinSize()) )
        {
          UInt uiRawBits = (2 * g_bitDepthY + g_bitDepthC) * rpcBestCU->getWidth(0) * rpcBestCU->getHeight(0) / 2;
          UInt uiBestBits = rpcBestCU->getTotalBits();
          if((uiBestBits > uiRawBits) || (rpcBestCU->getTotalCost() > m_pcRdCost->calcRdCost(uiRawBits, 0)))
          {
            xCheckIntraPCM (rpcBestCU, rpcTempCU);
            rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );
          }
        }
        if (bIsLosslessMode)
        {
          iQP = iMinQP;
        }
      }
    }

    m_pcEntropyCoder->resetBits();
    m_pcEntropyCoder->encodeSplitFlag( rpcBestCU, 0, uiDepth, true );
    rpcBestCU->getTotalBits() += m_pcEntropyCoder->getNumberOfWrittenBits(); // split bits
    rpcBestCU->getTotalBins() += ((TEncBinCABAC *)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
    rpcBestCU->getTotalCost()  = m_pcRdCost->calcRdCost( rpcBestCU->getTotalBits(), rpcBestCU->getTotalDistortion() );

    // Early CU determination
    if( m_pcEncCfg->getUseEarlyCU() && rpcBestCU->isSkipped(0) )
    {
      bSubBranch = false;
    }
    else
    {
      bSubBranch = true;
    }
  }
  else if(!(bSliceEnd && bInsidePicture))
  {
    bBoundary = true;
  }

  // copy orginal YUV samples to PCM buffer
  if( rpcBestCU->isLosslessCoded(0) && (rpcBestCU->getIPCMFlag(0) == false))
  {
    xFillPCMBuffer(rpcBestCU, m_ppcOrigYuv[uiDepth]);
  }
  if( (g_uiMaxCUWidth>>uiDepth) == rpcTempCU->getSlice()->getPPS()->getMinCuDQPSize() )
  {
    Int idQP = m_pcEncCfg->getMaxDeltaQP();
    iMinQP = Clip3( -rpcTempCU->getSlice()->getSPS()->getQpBDOffsetY(), MAX_QP, iBaseQP-idQP );
    iMaxQP = Clip3( -rpcTempCU->getSlice()->getSPS()->getQpBDOffsetY(), MAX_QP, iBaseQP+idQP );
  }
  else if( (g_uiMaxCUWidth>>uiDepth) > rpcTempCU->getSlice()->getPPS()->getMinCuDQPSize() )
  {
    iMinQP = iBaseQP;
    iMaxQP = iBaseQP;
  }
  else
  {
    Int iStartQP;
    if( pcPic->getCU( rpcTempCU->getAddr() )->getSliceSegmentStartCU(rpcTempCU->getZorderIdxInCU()) == pcSlice->getSliceSegmentCurStartCUAddr())
    {
      iStartQP = rpcTempCU->getQP(0);
    }
    else
    {
      UInt uiCurSliceStartPartIdx = pcSlice->getSliceSegmentCurStartCUAddr() % pcPic->getNumPartInCU() - rpcTempCU->getZorderIdxInCU();
      iStartQP = rpcTempCU->getQP(uiCurSliceStartPartIdx);
    }
    iMinQP = iStartQP;
    iMaxQP = iStartQP;
  }
  if ( m_pcEncCfg->getUseRateCtrl() )
  {
    iMinQP = m_pcRateCtrl->getRCQP();
    iMaxQP = m_pcRateCtrl->getRCQP();
  }

  if ( m_pcEncCfg->getCUTransquantBypassFlagForceValue() )
  {
    iMaxQP = iMinQP; // If all blocks are forced into using transquant bypass, do not loop here.
  }
    

  for (Int iQP=iMinQP; iQP<=iMaxQP; iQP++)
  {
    const Bool bIsLosslessMode = false; // False at this level. Next level down may set it to true.
    rpcTempCU->initEstData( uiDepth, iQP, bIsLosslessMode );

    // further split
    if( bSubBranch && uiDepth < g_uiMaxCUDepth - g_uiAddCUDepth )
    {
      UChar       uhNextDepth         = uiDepth+1;
      TComDataCU* pcSubBestPartCU     = m_ppcBestCU[uhNextDepth];
      TComDataCU* pcSubTempPartCU     = m_ppcTempCU[uhNextDepth];
        
        

      for ( UInt uiPartUnitIdx = 0; uiPartUnitIdx < 4; uiPartUnitIdx++ )
      {
        pcSubBestPartCU->initSubCU( rpcTempCU, uiPartUnitIdx, uhNextDepth, iQP );           // clear sub partition datas or init.
        pcSubTempPartCU->initSubCU( rpcTempCU, uiPartUnitIdx, uhNextDepth, iQP );           // clear sub partition datas or init.

        Bool bInSlice = pcSubBestPartCU->getSCUAddr()+pcSubBestPartCU->getTotalNumPart()>pcSlice->getSliceSegmentCurStartCUAddr()&&pcSubBestPartCU->getSCUAddr()<pcSlice->getSliceSegmentCurEndCUAddr();
        if(bInSlice && ( pcSubBestPartCU->getCUPelX() < pcSlice->getSPS()->getPicWidthInLumaSamples() ) && ( pcSubBestPartCU->getCUPelY() < pcSlice->getSPS()->getPicHeightInLumaSamples() ) )
        {
          if ( 0 == uiPartUnitIdx) //initialize RD with previous depth buffer
          {
            m_pppcRDSbacCoder[uhNextDepth][CI_CURR_BEST]->load(m_pppcRDSbacCoder[uiDepth][CI_CURR_BEST]);
          }
          else
          {
            m_pppcRDSbacCoder[uhNextDepth][CI_CURR_BEST]->load(m_pppcRDSbacCoder[uhNextDepth][CI_NEXT_BEST]);
          }

#if AMP_ENC_SPEEDUP
          if ( rpcBestCU->isIntra(0) )
          {
            xCompressCU( pcSubBestPartCU, pcSubTempPartCU, uhNextDepth, SIZE_NONE );
          }
          else
          {
            xCompressCU( pcSubBestPartCU, pcSubTempPartCU, uhNextDepth, rpcBestCU->getPartitionSize(0) );
          }
#else
            if ( rpcBestCU->isIntra(0) )
            {
                xCompressCU( pcSubBestPartCU, pcSubTempPartCU, uhNextDepth, SIZE_NONE );
            }
            else
            {
                xCompressCU( pcSubBestPartCU, pcSubTempPartCU, uhNextDepth, rpcBestCU->getPartitionSize(0) );
            }
            
            

            //xCompressCU( pcSubBestPartCU, pcSubTempPartCU, uhNextDepth );
#endif
            
          rpcTempCU->copyPartFrom( pcSubBestPartCU, uiPartUnitIdx, uhNextDepth );         // Keep best part data to current temporary data.
          xCopyYuv2Tmp( pcSubBestPartCU->getTotalNumPart()*uiPartUnitIdx, uhNextDepth );
        }
        else if (bInSlice)
        {
          pcSubBestPartCU->copyToPic( uhNextDepth );
          rpcTempCU->copyPartFrom( pcSubBestPartCU, uiPartUnitIdx, uhNextDepth );
        }
      } // end of for loop

      if( !bBoundary )
      {
        m_pcEntropyCoder->resetBits();
        m_pcEntropyCoder->encodeSplitFlag( rpcTempCU, 0, uiDepth, true );

        rpcTempCU->getTotalBits() += m_pcEntropyCoder->getNumberOfWrittenBits(); // split bits
        rpcTempCU->getTotalBins() += ((TEncBinCABAC *)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
      }
      rpcTempCU->getTotalCost()  = m_pcRdCost->calcRdCost( rpcTempCU->getTotalBits(), rpcTempCU->getTotalDistortion() );

      if( (g_uiMaxCUWidth>>uiDepth) == rpcTempCU->getSlice()->getPPS()->getMinCuDQPSize() && rpcTempCU->getSlice()->getPPS()->getUseDQP())
      {
        Bool hasResidual = false;
        for( UInt uiBlkIdx = 0; uiBlkIdx < rpcTempCU->getTotalNumPart(); uiBlkIdx ++)
        {
          if( ( pcPic->getCU( rpcTempCU->getAddr() )->getSliceSegmentStartCU(uiBlkIdx+rpcTempCU->getZorderIdxInCU()) == rpcTempCU->getSlice()->getSliceSegmentCurStartCUAddr() ) && 
              ( rpcTempCU->getCbf( uiBlkIdx, TEXT_LUMA ) || rpcTempCU->getCbf( uiBlkIdx, TEXT_CHROMA_U ) || rpcTempCU->getCbf( uiBlkIdx, TEXT_CHROMA_V ) ) )
          {
            hasResidual = true;
            break;
          }
        }

        UInt uiTargetPartIdx;
        if ( pcPic->getCU( rpcTempCU->getAddr() )->getSliceSegmentStartCU(rpcTempCU->getZorderIdxInCU()) != pcSlice->getSliceSegmentCurStartCUAddr() )
        {
          uiTargetPartIdx = pcSlice->getSliceSegmentCurStartCUAddr() % pcPic->getNumPartInCU() - rpcTempCU->getZorderIdxInCU();
        }
        else
        {
          uiTargetPartIdx = 0;
        }
        if ( hasResidual )
        {
#if !RDO_WITHOUT_DQP_BITS
          m_pcEntropyCoder->resetBits();
          m_pcEntropyCoder->encodeQP( rpcTempCU, uiTargetPartIdx, false );
          rpcTempCU->getTotalBits() += m_pcEntropyCoder->getNumberOfWrittenBits(); // dQP bits
          rpcTempCU->getTotalBins() += ((TEncBinCABAC *)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
          rpcTempCU->getTotalCost()  = m_pcRdCost->calcRdCost( rpcTempCU->getTotalBits(), rpcTempCU->getTotalDistortion() );
#endif

          Bool foundNonZeroCbf = false;
          rpcTempCU->setQPSubCUs( rpcTempCU->getRefQP( uiTargetPartIdx ), rpcTempCU, 0, uiDepth, foundNonZeroCbf );
          assert( foundNonZeroCbf );
        }
        else
        {
          rpcTempCU->setQPSubParts( rpcTempCU->getRefQP( uiTargetPartIdx ), 0, uiDepth ); // set QP to default QP
        }
      } // end of if

      m_pppcRDSbacCoder[uhNextDepth][CI_NEXT_BEST]->store(m_pppcRDSbacCoder[uiDepth][CI_TEMP_BEST]);

      Bool isEndOfSlice        = rpcBestCU->getSlice()->getSliceMode()==FIXED_NUMBER_OF_BYTES
                                 && (rpcBestCU->getTotalBits()>rpcBestCU->getSlice()->getSliceArgument()<<3);
      Bool isEndOfSliceSegment = rpcBestCU->getSlice()->getSliceSegmentMode()==FIXED_NUMBER_OF_BYTES
                                 && (rpcBestCU->getTotalBits()>rpcBestCU->getSlice()->getSliceSegmentArgument()<<3);
      if(isEndOfSlice||isEndOfSliceSegment)
      {
          
        rpcBestCU->getTotalCost()=rpcTempCU->getTotalCost()+1;
      }
        
     
      xCheckBestMode( rpcBestCU, rpcTempCU, uiDepth);                                               // RD compare current larger prediction
        
        //Counter 3
        int fnum =rpcBestCU->getSlice()->getPOC();
        if (fnum>1 && fnum<=input1 && numCTU<=91) //if(fnum>1 && fnum<=input1)
        {
             //printf("Checked in at Counter 3 \n\n");
            //printf("Depth = %d \n\n", rpcBestCU->getDepth(0));
            check1++;
            if(check1%21==0)
                
            {
                                    counter++;
                    //printf("counter = %d \n\n", counter);
                if(counter<=91)
                    xRoughMotionEsti(rpcBestCU, counter );
                
            }
        
        }
    }    // end of if                                                                              // with sub partitioned prediction.
      
  } // end of for loop QP
   
      rpcBestCU->copyToPic(uiDepth);                                                     // Copy Best data to Picture for next partition prediction.

  xCopyYuv2Pic( rpcBestCU->getPic(), rpcBestCU->getAddr(), rpcBestCU->getZorderIdxInCU(), uiDepth, uiDepth, rpcBestCU, uiLPelX, uiTPelY );   // Copy Yuv data to picture Yuv
  if( bBoundary ||(bSliceEnd && bInsidePicture))
  {
    return;
  }

  // Assert if Best prediction mode is NONE
  // Selected mode's RD-cost must be not MAX_DOUBLE.
    assert( rpcBestCU->getPartitionSize ( 0 ) != SIZE_NONE  );
    assert( rpcBestCU->getPredictionMode( 0 ) != MODE_NONE  );
    assert( rpcBestCU->getTotalCost     (   ) != MAX_DOUBLE );
    
}


/** finish encoding a cu and handle end-of-slice conditions
 * \param pcCU
 * \param uiAbsPartIdx
 * \param uiDepth 
 * \returns Void
 */
Void TEncCu::finishCU( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{
  TComPic* pcPic = pcCU->getPic();
  TComSlice * pcSlice = pcCU->getPic()->getSlice(pcCU->getPic()->getCurrSliceIdx());

  //Calculate end address
  UInt uiCUAddr = pcCU->getSCUAddr()+uiAbsPartIdx;

  UInt uiInternalAddress = pcPic->getPicSym()->getPicSCUAddr(pcSlice->getSliceSegmentCurEndCUAddr()-1) % pcPic->getNumPartInCU();
  UInt uiExternalAddress = pcPic->getPicSym()->getPicSCUAddr(pcSlice->getSliceSegmentCurEndCUAddr()-1) / pcPic->getNumPartInCU();
  UInt uiPosX = ( uiExternalAddress % pcPic->getFrameWidthInCU() ) * g_uiMaxCUWidth+ g_auiRasterToPelX[ g_auiZscanToRaster[uiInternalAddress] ];
  UInt uiPosY = ( uiExternalAddress / pcPic->getFrameWidthInCU() ) * g_uiMaxCUHeight+ g_auiRasterToPelY[ g_auiZscanToRaster[uiInternalAddress] ];
  UInt uiWidth = pcSlice->getSPS()->getPicWidthInLumaSamples();
  UInt uiHeight = pcSlice->getSPS()->getPicHeightInLumaSamples();
  while(uiPosX>=uiWidth||uiPosY>=uiHeight)
  {
    uiInternalAddress--;
    uiPosX = ( uiExternalAddress % pcPic->getFrameWidthInCU() ) * g_uiMaxCUWidth+ g_auiRasterToPelX[ g_auiZscanToRaster[uiInternalAddress] ];
    uiPosY = ( uiExternalAddress / pcPic->getFrameWidthInCU() ) * g_uiMaxCUHeight+ g_auiRasterToPelY[ g_auiZscanToRaster[uiInternalAddress] ];
  }
  uiInternalAddress++;
  if(uiInternalAddress==pcCU->getPic()->getNumPartInCU())
  {
    uiInternalAddress = 0;
    uiExternalAddress = pcPic->getPicSym()->getCUOrderMap(pcPic->getPicSym()->getInverseCUOrderMap(uiExternalAddress)+1);
  }
  UInt uiRealEndAddress = pcPic->getPicSym()->getPicSCUEncOrder(uiExternalAddress*pcPic->getNumPartInCU()+uiInternalAddress);

  // Encode slice finish
  Bool bTerminateSlice = false;
  if (uiCUAddr+(pcCU->getPic()->getNumPartInCU()>>(uiDepth<<1)) == uiRealEndAddress)
  {
    bTerminateSlice = true;
  }
  UInt uiGranularityWidth = g_uiMaxCUWidth;
  uiPosX = pcCU->getCUPelX() + g_auiRasterToPelX[ g_auiZscanToRaster[uiAbsPartIdx] ];
  uiPosY = pcCU->getCUPelY() + g_auiRasterToPelY[ g_auiZscanToRaster[uiAbsPartIdx] ];
  Bool granularityBoundary=((uiPosX+pcCU->getWidth(uiAbsPartIdx))%uiGranularityWidth==0||(uiPosX+pcCU->getWidth(uiAbsPartIdx)==uiWidth))
    &&((uiPosY+pcCU->getHeight(uiAbsPartIdx))%uiGranularityWidth==0||(uiPosY+pcCU->getHeight(uiAbsPartIdx)==uiHeight));
  
  if(granularityBoundary)
  {
    // The 1-terminating bit is added to all streams, so don't add it here when it's 1.
    if (!bTerminateSlice)
      m_pcEntropyCoder->encodeTerminatingBit( bTerminateSlice ? 1 : 0 );
  }
  
  Int numberOfWrittenBits = 0;
  if (m_pcBitCounter)
  {
    numberOfWrittenBits = m_pcEntropyCoder->getNumberOfWrittenBits();
  }
  
  // Calculate slice end IF this CU puts us over slice bit size.
  UInt iGranularitySize = pcCU->getPic()->getNumPartInCU();
  Int iGranularityEnd = ((pcCU->getSCUAddr()+uiAbsPartIdx)/iGranularitySize)*iGranularitySize;
  if(iGranularityEnd<=pcSlice->getSliceSegmentCurStartCUAddr()) 
  {
    iGranularityEnd+=max(iGranularitySize,(pcCU->getPic()->getNumPartInCU()>>(uiDepth<<1)));
  }
  // Set slice end parameter
  if(pcSlice->getSliceMode()==FIXED_NUMBER_OF_BYTES&&!pcSlice->getFinalized()&&pcSlice->getSliceBits()+numberOfWrittenBits>pcSlice->getSliceArgument()<<3) 
  {
    pcSlice->setSliceSegmentCurEndCUAddr(iGranularityEnd);
    pcSlice->setSliceCurEndCUAddr(iGranularityEnd);
    return;
  }
  // Set dependent slice end parameter
  if(pcSlice->getSliceSegmentMode()==FIXED_NUMBER_OF_BYTES&&!pcSlice->getFinalized()&&pcSlice->getSliceSegmentBits()+numberOfWrittenBits > pcSlice->getSliceSegmentArgument()<<3) 
  {
    pcSlice->setSliceSegmentCurEndCUAddr(iGranularityEnd);
    return;
  }
  if(granularityBoundary)
  {
    pcSlice->setSliceBits( (UInt)(pcSlice->getSliceBits() + numberOfWrittenBits) );
    pcSlice->setSliceSegmentBits(pcSlice->getSliceSegmentBits()+numberOfWrittenBits);
    if (m_pcBitCounter)
    {
      m_pcEntropyCoder->resetBits();      
    }
  }
}

/** Compute QP for each CU
 * \param pcCU Target CU
 * \param uiDepth CU depth
 * \returns quantization parameter
 */
Int TEncCu::xComputeQP( TComDataCU* pcCU, UInt uiDepth )
{
  Int iBaseQp = pcCU->getSlice()->getSliceQp();
  Int iQpOffset = 0;
  if ( m_pcEncCfg->getUseAdaptiveQP() )
  {
    TEncPic* pcEPic = dynamic_cast<TEncPic*>( pcCU->getPic() );
    UInt uiAQDepth = min( uiDepth, pcEPic->getMaxAQDepth()-1 );
    TEncPicQPAdaptationLayer* pcAQLayer = pcEPic->getAQLayer( uiAQDepth );
    UInt uiAQUPosX = pcCU->getCUPelX() / pcAQLayer->getAQPartWidth();
    UInt uiAQUPosY = pcCU->getCUPelY() / pcAQLayer->getAQPartHeight();
    UInt uiAQUStride = pcAQLayer->getAQPartStride();
    TEncQPAdaptationUnit* acAQU = pcAQLayer->getQPAdaptationUnit();

    Double dMaxQScale = pow(2.0, m_pcEncCfg->getQPAdaptationRange()/6.0);
    Double dAvgAct = pcAQLayer->getAvgActivity();
    Double dCUAct = acAQU[uiAQUPosY * uiAQUStride + uiAQUPosX].getActivity();
    Double dNormAct = (dMaxQScale*dCUAct + dAvgAct) / (dCUAct + dMaxQScale*dAvgAct);
    Double dQpOffset = log(dNormAct) / log(2.0) * 6.0;
    iQpOffset = Int(floor( dQpOffset + 0.49999 ));
  }
  return Clip3(-pcCU->getSlice()->getSPS()->getQpBDOffsetY(), MAX_QP, iBaseQp+iQpOffset );
}

/** encode a CU block recursively
 * \param pcCU
 * \param uiAbsPartIdx
 * \param uiDepth 
 * \returns Void
 */
Void TEncCu::xEncodeCU( TComDataCU* pcCU, UInt uiAbsPartIdx, UInt uiDepth )
{
  TComPic* pcPic = pcCU->getPic();
  
  Bool bBoundary = false;
  UInt uiLPelX   = pcCU->getCUPelX() + g_auiRasterToPelX[ g_auiZscanToRaster[uiAbsPartIdx] ];
  UInt uiRPelX   = uiLPelX + (g_uiMaxCUWidth>>uiDepth)  - 1;
  UInt uiTPelY   = pcCU->getCUPelY() + g_auiRasterToPelY[ g_auiZscanToRaster[uiAbsPartIdx] ];
  UInt uiBPelY   = uiTPelY + (g_uiMaxCUHeight>>uiDepth) - 1;
  
  TComSlice * pcSlice = pcCU->getPic()->getSlice(pcCU->getPic()->getCurrSliceIdx());
  // If slice start is within this cu...
  Bool bSliceStart = pcSlice->getSliceSegmentCurStartCUAddr() > pcPic->getPicSym()->getInverseCUOrderMap(pcCU->getAddr())*pcCU->getPic()->getNumPartInCU()+uiAbsPartIdx && 
    pcSlice->getSliceSegmentCurStartCUAddr() < pcPic->getPicSym()->getInverseCUOrderMap(pcCU->getAddr())*pcCU->getPic()->getNumPartInCU()+uiAbsPartIdx+( pcPic->getNumPartInCU() >> (uiDepth<<1) );
  // We need to split, so don't try these modes.
  if(!bSliceStart&&( uiRPelX < pcSlice->getSPS()->getPicWidthInLumaSamples() ) && ( uiBPelY < pcSlice->getSPS()->getPicHeightInLumaSamples() ) )
  {
    m_pcEntropyCoder->encodeSplitFlag( pcCU, uiAbsPartIdx, uiDepth );
  }
  else
  {
    bBoundary = true;
  }
  
  if( ( ( uiDepth < pcCU->getDepth( uiAbsPartIdx ) ) && ( uiDepth < (g_uiMaxCUDepth-g_uiAddCUDepth) ) ) || bBoundary )
  {
    UInt uiQNumParts = ( pcPic->getNumPartInCU() >> (uiDepth<<1) )>>2;
    if( (g_uiMaxCUWidth>>uiDepth) == pcCU->getSlice()->getPPS()->getMinCuDQPSize() && pcCU->getSlice()->getPPS()->getUseDQP())
    {
      setdQPFlag(true);
    }
    for ( UInt uiPartUnitIdx = 0; uiPartUnitIdx < 4; uiPartUnitIdx++, uiAbsPartIdx+=uiQNumParts )
    {
      uiLPelX   = pcCU->getCUPelX() + g_auiRasterToPelX[ g_auiZscanToRaster[uiAbsPartIdx] ];
      uiTPelY   = pcCU->getCUPelY() + g_auiRasterToPelY[ g_auiZscanToRaster[uiAbsPartIdx] ];
      Bool bInSlice = pcCU->getSCUAddr()+uiAbsPartIdx+uiQNumParts>pcSlice->getSliceSegmentCurStartCUAddr()&&pcCU->getSCUAddr()+uiAbsPartIdx<pcSlice->getSliceSegmentCurEndCUAddr();
      if(bInSlice&&( uiLPelX < pcSlice->getSPS()->getPicWidthInLumaSamples() ) && ( uiTPelY < pcSlice->getSPS()->getPicHeightInLumaSamples() ) )
      {
        xEncodeCU( pcCU, uiAbsPartIdx, uiDepth+1 );
      }
    }
    return;
  }
  
  if( (g_uiMaxCUWidth>>uiDepth) >= pcCU->getSlice()->getPPS()->getMinCuDQPSize() && pcCU->getSlice()->getPPS()->getUseDQP())
  {
    setdQPFlag(true);
  }
  if (pcCU->getSlice()->getPPS()->getTransquantBypassEnableFlag())
  {
    m_pcEntropyCoder->encodeCUTransquantBypassFlag( pcCU, uiAbsPartIdx );
  }
  if( !pcCU->getSlice()->isIntra() )
  {
    m_pcEntropyCoder->encodeSkipFlag( pcCU, uiAbsPartIdx );
  }
  
  if( pcCU->isSkipped( uiAbsPartIdx ) )
  {
    m_pcEntropyCoder->encodeMergeIndex( pcCU, uiAbsPartIdx );
    finishCU(pcCU,uiAbsPartIdx,uiDepth);
    return;
  }
  m_pcEntropyCoder->encodePredMode( pcCU, uiAbsPartIdx );
  
  m_pcEntropyCoder->encodePartSize( pcCU, uiAbsPartIdx, uiDepth );
  
  if (pcCU->isIntra( uiAbsPartIdx ) && pcCU->getPartitionSize( uiAbsPartIdx ) == SIZE_2Nx2N )
  {
    m_pcEntropyCoder->encodeIPCMInfo( pcCU, uiAbsPartIdx );

    if(pcCU->getIPCMFlag(uiAbsPartIdx))
    {
      // Encode slice finish
      finishCU(pcCU,uiAbsPartIdx,uiDepth);
      return;
    }
  }

  // prediction Info ( Intra : direction mode, Inter : Mv, reference idx )
  m_pcEntropyCoder->encodePredInfo( pcCU, uiAbsPartIdx );
  
  // Encode Coefficients
  Bool bCodeDQP = getdQPFlag();
  m_pcEntropyCoder->encodeCoeff( pcCU, uiAbsPartIdx, uiDepth, pcCU->getWidth (uiAbsPartIdx), pcCU->getHeight(uiAbsPartIdx), bCodeDQP );
  setdQPFlag( bCodeDQP );

  // --- write terminating bit ---
  finishCU(pcCU,uiAbsPartIdx,uiDepth);
}

Int xCalcHADs8x8_ISlice(Pel *piOrg, Int iStrideOrg) 
{
  Int k, i, j, jj;
  Int diff[64], m1[8][8], m2[8][8], m3[8][8], iSumHad = 0;

  for( k = 0; k < 64; k += 8 )
  {
    diff[k+0] = piOrg[0] ;
    diff[k+1] = piOrg[1] ;
    diff[k+2] = piOrg[2] ;
    diff[k+3] = piOrg[3] ;
    diff[k+4] = piOrg[4] ;
    diff[k+5] = piOrg[5] ;
    diff[k+6] = piOrg[6] ;
    diff[k+7] = piOrg[7] ;
 
    piOrg += iStrideOrg;
  }
  
  //horizontal
  for (j=0; j < 8; j++)
  {
    jj = j << 3;
    m2[j][0] = diff[jj  ] + diff[jj+4];
    m2[j][1] = diff[jj+1] + diff[jj+5];
    m2[j][2] = diff[jj+2] + diff[jj+6];
    m2[j][3] = diff[jj+3] + diff[jj+7];
    m2[j][4] = diff[jj  ] - diff[jj+4];
    m2[j][5] = diff[jj+1] - diff[jj+5];
    m2[j][6] = diff[jj+2] - diff[jj+6];
    m2[j][7] = diff[jj+3] - diff[jj+7];
    
    m1[j][0] = m2[j][0] + m2[j][2];
    m1[j][1] = m2[j][1] + m2[j][3];
    m1[j][2] = m2[j][0] - m2[j][2];
    m1[j][3] = m2[j][1] - m2[j][3];
    m1[j][4] = m2[j][4] + m2[j][6];
    m1[j][5] = m2[j][5] + m2[j][7];
    m1[j][6] = m2[j][4] - m2[j][6];
    m1[j][7] = m2[j][5] - m2[j][7];
    
    m2[j][0] = m1[j][0] + m1[j][1];
    m2[j][1] = m1[j][0] - m1[j][1];
    m2[j][2] = m1[j][2] + m1[j][3];
    m2[j][3] = m1[j][2] - m1[j][3];
    m2[j][4] = m1[j][4] + m1[j][5];
    m2[j][5] = m1[j][4] - m1[j][5];
    m2[j][6] = m1[j][6] + m1[j][7];
    m2[j][7] = m1[j][6] - m1[j][7];
  }
  
  //vertical
  for (i=0; i < 8; i++)
  {
    m3[0][i] = m2[0][i] + m2[4][i];
    m3[1][i] = m2[1][i] + m2[5][i];
    m3[2][i] = m2[2][i] + m2[6][i];
    m3[3][i] = m2[3][i] + m2[7][i];
    m3[4][i] = m2[0][i] - m2[4][i];
    m3[5][i] = m2[1][i] - m2[5][i];
    m3[6][i] = m2[2][i] - m2[6][i];
    m3[7][i] = m2[3][i] - m2[7][i];
    
    m1[0][i] = m3[0][i] + m3[2][i];
    m1[1][i] = m3[1][i] + m3[3][i];
    m1[2][i] = m3[0][i] - m3[2][i];
    m1[3][i] = m3[1][i] - m3[3][i];
    m1[4][i] = m3[4][i] + m3[6][i];
    m1[5][i] = m3[5][i] + m3[7][i];
    m1[6][i] = m3[4][i] - m3[6][i];
    m1[7][i] = m3[5][i] - m3[7][i];
    
    m2[0][i] = m1[0][i] + m1[1][i];
    m2[1][i] = m1[0][i] - m1[1][i];
    m2[2][i] = m1[2][i] + m1[3][i];
    m2[3][i] = m1[2][i] - m1[3][i];
    m2[4][i] = m1[4][i] + m1[5][i];
    m2[5][i] = m1[4][i] - m1[5][i];
    m2[6][i] = m1[6][i] + m1[7][i];
    m2[7][i] = m1[6][i] - m1[7][i];
  }
  
  for (i = 0; i < 8; i++)
  {
    for (j = 0; j < 8; j++)
    {
      iSumHad += abs(m2[i][j]);
    }
  }
  iSumHad -= abs(m2[0][0]);
  iSumHad =(iSumHad+2)>>2;
  return(iSumHad);
}

Int  TEncCu::updateLCUDataISlice(TComDataCU* pcCU, Int LCUIdx, Int width, Int height)
{
  Int  xBl, yBl; 
  const Int iBlkSize = 8;

  Pel* pOrgInit   = pcCU->getPic()->getPicYuvOrg()->getLumaAddr(pcCU->getAddr(), 0);
  Int  iStrideOrig = pcCU->getPic()->getPicYuvOrg()->getStride();
  Pel  *pOrg;

  Int iSumHad = 0;
  for ( yBl=0; (yBl+iBlkSize)<=height; yBl+= iBlkSize)
  {
    for ( xBl=0; (xBl+iBlkSize)<=width; xBl+= iBlkSize)
    {
      pOrg = pOrgInit + iStrideOrig*yBl + xBl; 
      iSumHad += xCalcHADs8x8_ISlice(pOrg, iStrideOrig);
    }
  }
  return(iSumHad);
}

/** check RD costs for a CU block encoded with merge
 * \param rpcBestCU
 * \param rpcTempCU
 * \returns Void
 */

Void TEncCu::xCheckRDCostMerge2Nx2N( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, Bool *earlyDetectionSkipMode )
{
  assert( rpcTempCU->getSlice()->getSliceType() != I_SLICE );
  TComMvField  cMvFieldNeighbours[ 2 * MRG_MAX_NUM_CANDS ]; // double length for mv of both lists
  UChar uhInterDirNeighbours[MRG_MAX_NUM_CANDS];
  Int numValidMergeCand = 0;
  const Bool bTransquantBypassFlag = rpcTempCU->getCUTransquantBypass(0);

  for( UInt ui = 0; ui < rpcTempCU->getSlice()->getMaxNumMergeCand(); ++ui )
  {
    uhInterDirNeighbours[ui] = 0;
  }
  UChar uhDepth = rpcTempCU->getDepth( 0 );
  rpcTempCU->setPartSizeSubParts( SIZE_2Nx2N, 0, uhDepth ); // interprets depth relative to LCU level
  rpcTempCU->getInterMergeCandidates( 0, 0, cMvFieldNeighbours,uhInterDirNeighbours, numValidMergeCand );
  
  Int mergeCandBuffer[MRG_MAX_NUM_CANDS];
  for( UInt ui = 0; ui < numValidMergeCand; ++ui )
  {
    mergeCandBuffer[ui] = 0;
  }

  Bool bestIsSkip = false;

  UInt iteration;
  if ( rpcTempCU->isLosslessCoded(0))
  {
    iteration = 1;
  }
  else 
  {
    iteration = 2;
  }

  for( UInt uiNoResidual = 0; uiNoResidual < iteration; ++uiNoResidual )
  {
    for( UInt uiMergeCand = 0; uiMergeCand < numValidMergeCand; ++uiMergeCand )
    {
      if(!(uiNoResidual==1 && mergeCandBuffer[uiMergeCand]==1))
      {
        if( !(bestIsSkip && uiNoResidual == 0) )
        {
          // set MC parameters
          rpcTempCU->setPredModeSubParts( MODE_INTER, 0, uhDepth ); // interprets depth relative to LCU level
          rpcTempCU->setCUTransquantBypassSubParts( bTransquantBypassFlag,     0, uhDepth );
          rpcTempCU->setPartSizeSubParts( SIZE_2Nx2N, 0, uhDepth ); // interprets depth relative to LCU level
          rpcTempCU->setMergeFlagSubParts( true, 0, 0, uhDepth ); // interprets depth relative to LCU level
          rpcTempCU->setMergeIndexSubParts( uiMergeCand, 0, 0, uhDepth ); // interprets depth relative to LCU level
          rpcTempCU->setInterDirSubParts( uhInterDirNeighbours[uiMergeCand], 0, 0, uhDepth ); // interprets depth relative to LCU level
          rpcTempCU->getCUMvField( REF_PIC_LIST_0 )->setAllMvField( cMvFieldNeighbours[0 + 2*uiMergeCand], SIZE_2Nx2N, 0, 0 ); // interprets depth relative to rpcTempCU level
          rpcTempCU->getCUMvField( REF_PIC_LIST_1 )->setAllMvField( cMvFieldNeighbours[1 + 2*uiMergeCand], SIZE_2Nx2N, 0, 0 ); // interprets depth relative to rpcTempCU level
          
          // do MC
          m_pcPredSearch->motionCompensation ( rpcTempCU, m_ppcPredYuvTemp[uhDepth] );
          // estimate residual and encode everything
          m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU,
                                                    m_ppcOrigYuv    [uhDepth],
                                                    m_ppcPredYuvTemp[uhDepth],
                                                    m_ppcResiYuvTemp[uhDepth],
                                                    m_ppcResiYuvBest[uhDepth],
                                                    m_ppcRecoYuvTemp[uhDepth],
                                                    (uiNoResidual? true:false));
          
          
          if ( uiNoResidual == 0 && rpcTempCU->getQtRootCbf(0) == 0 )
          {
            // If no residual when allowing for one, then set mark to not try case where residual is forced to 0
            mergeCandBuffer[uiMergeCand] = 1;
          }
          
          rpcTempCU->setSkipFlagSubParts( rpcTempCU->getQtRootCbf(0) == 0, 0, uhDepth );
          Int orgQP = rpcTempCU->getQP( 0 );
          xCheckDQP( rpcTempCU );
          xCheckBestMode(rpcBestCU, rpcTempCU, uhDepth);
          rpcTempCU->initEstData( uhDepth, orgQP, bTransquantBypassFlag );
          
          if( m_pcEncCfg->getUseFastDecisionForMerge() && !bestIsSkip )
          {
            bestIsSkip = rpcBestCU->getQtRootCbf(0) == 0;
          }
        }
      }
    }
    
    if(uiNoResidual == 0 && m_pcEncCfg->getUseEarlySkipDetection())
    {
      if(rpcBestCU->getQtRootCbf( 0 ) == 0)
      {
        if( rpcBestCU->getMergeFlag( 0 ))
        {
          *earlyDetectionSkipMode = true;
        }
        else
        {
          Int absoulte_MV=0;
          for ( UInt uiRefListIdx = 0; uiRefListIdx < 2; uiRefListIdx++ )
          {
            if ( rpcBestCU->getSlice()->getNumRefIdx( RefPicList( uiRefListIdx ) ) > 0 )
            {
              TComCUMvField* pcCUMvField = rpcBestCU->getCUMvField(RefPicList( uiRefListIdx ));
              Int iHor = pcCUMvField->getMvd( 0 ).getAbsHor();
              Int iVer = pcCUMvField->getMvd( 0 ).getAbsVer();
              absoulte_MV+=iHor+iVer;
            }
          }
          
          if(absoulte_MV == 0)
          {
            *earlyDetectionSkipMode = true;
          }
        }
      }
    }
  }
}


#if AMP_MRG
Void TEncCu::xCheckRDCostInter( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, PartSize ePartSize, Bool bUseMRG)
#else
Void TEncCu::xCheckRDCostInter( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, PartSize ePartSize, int uiDepth, int signal1 )
#endif
{
  UChar uhDepth = rpcTempCU->getDepth( 0 );
  
  rpcTempCU->setDepthSubParts( uhDepth, 0 );
  
  rpcTempCU->setSkipFlagSubParts( false, 0, uhDepth );

  rpcTempCU->setPartSizeSubParts  ( ePartSize,  0, uhDepth );
  rpcTempCU->setPredModeSubParts  ( MODE_INTER, 0, uhDepth );
    
    
#if AMP_MRG
  rpcTempCU->setMergeAMP (true);
  m_pcPredSearch->predInterSearch ( rpcTempCU, m_ppcOrigYuv[uhDepth], m_ppcPredYuvTemp[uhDepth], m_ppcResiYuvTemp[uhDepth], m_ppcRecoYuvTemp[uhDepth], false, bUseMRG );
#else
  m_pcPredSearch->predInterSearch ( rpcTempCU, m_ppcOrigYuv[uhDepth], m_ppcPredYuvTemp[uhDepth], m_ppcResiYuvTemp[uhDepth], m_ppcRecoYuvTemp[uhDepth] );
#endif

#if AMP_MRG
  if ( !rpcTempCU->getMergeAMP() )
  {
    return;
  }
#endif

  m_pcPredSearch->encodeResAndCalcRdInterCU( rpcTempCU, m_ppcOrigYuv[uhDepth], m_ppcPredYuvTemp[uhDepth], m_ppcResiYuvTemp[uhDepth], m_ppcResiYuvBest[uhDepth], m_ppcRecoYuvTemp[uhDepth], false );
  rpcTempCU->getTotalCost()  = m_pcRdCost->calcRdCost( rpcTempCU->getTotalBits(), rpcTempCU->getTotalDistortion() );

  xCheckDQP( rpcTempCU );
  xCheckBestMode(rpcBestCU, rpcTempCU, uhDepth);
    
     // Rough Motion Estimation 2N x 2N (N = 32) (Y)
    // Counter 2
    //printf("Checked in at Counter 2 \n\n");
    
if(rpcBestCU->getSlice()->getPOC()>1 && rpcBestCU->getSlice()->getPOC()<=input1 && uiDepth==0 && signal1==1 && rpcBestCU->getSlice()->getSliceType() != I_SLICE && numCTU<=91) // select frame index and CTU size
    {
        Pel* Pred = m_ppcPredYuvTemp[uhDepth]->getLumaAddr();
        Pel* Org  = m_ppcOrigYuv[uhDepth]->getLumaAddr();
        flag2=0;
        for(int r=(d1)*64; r<64+(d1)*64; r++)
        {
            for(int c=(c1)*64; c<64+(c1)*64; c++)
            {
                flag2++;
                residualy[r][c]=abs(*Pred - *Org );
                residualy_new[r][c]=abs(*Pred - *Org );
                if(!(flag2==64*64))
                {
                    Pred++;
                    Org++;
                }
            }
        }
        ::numCTU++; // count number of CTUs for a given frame
            c1++;
            if(c1==13)
            {
                d1++;
                c1=0;
            }  //change row since only 13 CTUs can be accommodated here.
        
    }
  
    
  
}
Void TEncCu::xRoughMotionEsti( TComDataCU*& rpcBestCU,int vari )
{

    // Counter 4
    ::numCTU = vari;
    printf("NUMCTU %d \n\n", ::numCTU);
    int row = floor((::numCTU-1)/13); // 0  to 6
    int col = (::numCTU-1)%13; // 0 to 12
    int dep = rpcBestCU->getDepth(0);
    int bsize = pow(2, 6-dep);
    int csize = bsize/2;
    
    for (int i=0+row*64; i<64+row*64; i++)
    {
        for(int j=0+col*64; j<64+col*64;j++)
        {
            for(int integer=1; integer<200; integer++)
            {
                if(i==bsize*integer-1 || j==bsize*integer-1)
                residualy[i][j]=255;
                
                if(rpcBestCU->getPartitionSize(0)==SIZE_Nx2N)
                {
                    if(j==csize*integer-1)
                        residualy[i][j]=255;
                }
                else if(rpcBestCU->getPartitionSize(0)==SIZE_2NxN)
                {
                   if(i==csize*integer-1 )
                       residualy[i][j]=255;
                }
                else if(rpcBestCU->getPartitionSize(0)==SIZE_NxN)
                {
                    if(i==csize*integer-1 || j==csize*integer-1)
                       residualy[i][j]=255;
                }
                
            }
            
        }
    }
    
    /*if(::numCTU==91)
    {
        
        ofstream arquivoOut0("rd.dat");
        for(int t=1; t<=91; t++)
        {
            arquivoOut0<< ctucost[t]/10000<<"\t";
        }
        arquivoOut0.close();
    }*/
    
    if (::numCTU==91)
    {
        framecount++;
        sprintf (buffer, "%d best.dat", framecount);
        ofstream arquivoOut0(buffer);
        //ofstream arquivoOut0("y_best.dat");
        
        for (int i=0; i<64*7; i++)
        {
            for(int j=0; j<64*13; j++)
            {
                
                if(j!= 64*13-1) arquivoOut0<< residualy[i][j] << "\t";
                else arquivoOut0<< residualy[i][j] << "\n";
                
            }
        }
        arquivoOut0.close();
    }
    
    if (::numCTU==91)
    {
        sprintf (buffer1, "%d rough.dat", framecount);
        ofstream arquivoOut0(buffer1);
        
        //ofstream arquivoOut0("y_rough.dat");

        for (int i=0; i<64*7; i++)
        {
            for(int j=0; j<64*13; j++)
            {
                
                if(j!= 64*13-1) arquivoOut0<< residualy_new[i][j] << "\t";
                else arquivoOut0<< residualy_new[i][j] << "\n";
                
            }
        }
        arquivoOut0.close();
        for (int i=0; i<481; i++)
        {
            for (int j=0; j<833; j++)
            {
                residualy[i][j] = 0;
                residualy_new[i][j] = 0;
            }
        }
        
    }

    


}

Void TEncCu::xCheckRDCostIntra( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, PartSize eSize )
{
  UInt uiDepth = rpcTempCU->getDepth( 0 );
  
  rpcTempCU->setSkipFlagSubParts( false, 0, uiDepth );

  rpcTempCU->setPartSizeSubParts( eSize, 0, uiDepth );
  rpcTempCU->setPredModeSubParts( MODE_INTRA, 0, uiDepth );
  
  Bool bSeparateLumaChroma = true; // choose estimation mode
  UInt uiPreCalcDistC      = 0;
  if( !bSeparateLumaChroma )
  {
    m_pcPredSearch->preestChromaPredMode( rpcTempCU, m_ppcOrigYuv[uiDepth], m_ppcPredYuvTemp[uiDepth] );
  }
  m_pcPredSearch  ->estIntraPredQT      ( rpcTempCU, m_ppcOrigYuv[uiDepth], m_ppcPredYuvTemp[uiDepth], m_ppcResiYuvTemp[uiDepth], m_ppcRecoYuvTemp[uiDepth], uiPreCalcDistC, bSeparateLumaChroma );

  m_ppcRecoYuvTemp[uiDepth]->copyToPicLuma(rpcTempCU->getPic()->getPicYuvRec(), rpcTempCU->getAddr(), rpcTempCU->getZorderIdxInCU() );
  
  m_pcPredSearch  ->estIntraPredChromaQT( rpcTempCU, m_ppcOrigYuv[uiDepth], m_ppcPredYuvTemp[uiDepth], m_ppcResiYuvTemp[uiDepth], m_ppcRecoYuvTemp[uiDepth], uiPreCalcDistC );
  
  m_pcEntropyCoder->resetBits();
  if ( rpcTempCU->getSlice()->getPPS()->getTransquantBypassEnableFlag())
  {
    m_pcEntropyCoder->encodeCUTransquantBypassFlag( rpcTempCU, 0,          true );
  }
  m_pcEntropyCoder->encodeSkipFlag ( rpcTempCU, 0,          true );
  m_pcEntropyCoder->encodePredMode( rpcTempCU, 0,          true );
  m_pcEntropyCoder->encodePartSize( rpcTempCU, 0, uiDepth, true );
  m_pcEntropyCoder->encodePredInfo( rpcTempCU, 0,          true );
  m_pcEntropyCoder->encodeIPCMInfo(rpcTempCU, 0, true );

  // Encode Coefficients
  Bool bCodeDQP = getdQPFlag();
  m_pcEntropyCoder->encodeCoeff( rpcTempCU, 0, uiDepth, rpcTempCU->getWidth (0), rpcTempCU->getHeight(0), bCodeDQP );
  setdQPFlag( bCodeDQP );
  
  m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[uiDepth][CI_TEMP_BEST]);
  
  rpcTempCU->getTotalBits() = m_pcEntropyCoder->getNumberOfWrittenBits();
  rpcTempCU->getTotalBins() = ((TEncBinCABAC *)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
  rpcTempCU->getTotalCost() = m_pcRdCost->calcRdCost( rpcTempCU->getTotalBits(), rpcTempCU->getTotalDistortion() );
  
  xCheckDQP( rpcTempCU );
  xCheckBestMode(rpcBestCU, rpcTempCU, uiDepth);
}

/** Check R-D costs for a CU with PCM mode. 
 * \param rpcBestCU pointer to best mode CU data structure
 * \param rpcTempCU pointer to testing mode CU data structure
 * \returns Void
 * 
 * \note Current PCM implementation encodes sample values in a lossless way. The distortion of PCM mode CUs are zero. PCM mode is selected if the best mode yields bits greater than that of PCM mode.
 */
Void TEncCu::xCheckIntraPCM( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU )
{
  UInt uiDepth = rpcTempCU->getDepth( 0 );

  rpcTempCU->setSkipFlagSubParts( false, 0, uiDepth );

  rpcTempCU->setIPCMFlag(0, true);
  rpcTempCU->setIPCMFlagSubParts (true, 0, rpcTempCU->getDepth(0));
  rpcTempCU->setPartSizeSubParts( SIZE_2Nx2N, 0, uiDepth );
  rpcTempCU->setPredModeSubParts( MODE_INTRA, 0, uiDepth );
  rpcTempCU->setTrIdxSubParts ( 0, 0, uiDepth );

  m_pcPredSearch->IPCMSearch( rpcTempCU, m_ppcOrigYuv[uiDepth], m_ppcPredYuvTemp[uiDepth], m_ppcResiYuvTemp[uiDepth], m_ppcRecoYuvTemp[uiDepth]);

  m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiDepth][CI_CURR_BEST]);

  m_pcEntropyCoder->resetBits();
  if ( rpcTempCU->getSlice()->getPPS()->getTransquantBypassEnableFlag())
  {
    m_pcEntropyCoder->encodeCUTransquantBypassFlag( rpcTempCU, 0,          true );
  }
  m_pcEntropyCoder->encodeSkipFlag ( rpcTempCU, 0,          true );
  m_pcEntropyCoder->encodePredMode ( rpcTempCU, 0,          true );
  m_pcEntropyCoder->encodePartSize ( rpcTempCU, 0, uiDepth, true );
  m_pcEntropyCoder->encodeIPCMInfo ( rpcTempCU, 0, true );

  m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[uiDepth][CI_TEMP_BEST]);

  rpcTempCU->getTotalBits() = m_pcEntropyCoder->getNumberOfWrittenBits();
  rpcTempCU->getTotalBins() = ((TEncBinCABAC *)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
  rpcTempCU->getTotalCost() = m_pcRdCost->calcRdCost( rpcTempCU->getTotalBits(), rpcTempCU->getTotalDistortion() );

  xCheckDQP( rpcTempCU );
  xCheckBestMode( rpcBestCU, rpcTempCU, uiDepth);
    
    
}

/** check whether current try is the best with identifying the depth of current try
 * \param rpcBestCU
 * \param rpcTempCU
 * \returns Void
 */
Void TEncCu::xCheckBestMode( TComDataCU*& rpcBestCU, TComDataCU*& rpcTempCU, UInt uiDepth )
{
    if( rpcTempCU->getTotalCost() < rpcBestCU->getTotalCost() )
  {
    TComYuv* pcYuv;
    // Change Information data
    TComDataCU* pcCU = rpcBestCU;
    rpcBestCU = rpcTempCU;
    rpcTempCU = pcCU;

    // Change Prediction data
    pcYuv = m_ppcPredYuvBest[uiDepth];
    m_ppcPredYuvBest[uiDepth] = m_ppcPredYuvTemp[uiDepth];
    m_ppcPredYuvTemp[uiDepth] = pcYuv;

    // Change Reconstruction data
    pcYuv = m_ppcRecoYuvBest[uiDepth];
    m_ppcRecoYuvBest[uiDepth] = m_ppcRecoYuvTemp[uiDepth];
    m_ppcRecoYuvTemp[uiDepth] = pcYuv;

    pcYuv = NULL;
    pcCU  = NULL;

    // store temp best CI for next CU coding
    m_pppcRDSbacCoder[uiDepth][CI_TEMP_BEST]->store(m_pppcRDSbacCoder[uiDepth][CI_NEXT_BEST]);
      
  }
}

Void TEncCu::xCheckDQP( TComDataCU* pcCU )
{
  UInt uiDepth = pcCU->getDepth( 0 );

  if( pcCU->getSlice()->getPPS()->getUseDQP() && (g_uiMaxCUWidth>>uiDepth) >= pcCU->getSlice()->getPPS()->getMinCuDQPSize() )
  {
    if ( pcCU->getCbf( 0, TEXT_LUMA, 0 ) || pcCU->getCbf( 0, TEXT_CHROMA_U, 0 ) || pcCU->getCbf( 0, TEXT_CHROMA_V, 0 ) )
    {
#if !RDO_WITHOUT_DQP_BITS
      m_pcEntropyCoder->resetBits();
      m_pcEntropyCoder->encodeQP( pcCU, 0, false );
      pcCU->getTotalBits() += m_pcEntropyCoder->getNumberOfWrittenBits(); // dQP bits
      pcCU->getTotalBins() += ((TEncBinCABAC *)((TEncSbac*)m_pcEntropyCoder->m_pcEntropyCoderIf)->getEncBinIf())->getBinsCoded();
      pcCU->getTotalCost() = m_pcRdCost->calcRdCost( pcCU->getTotalBits(), pcCU->getTotalDistortion() );
#endif
    }
    else
    {
      pcCU->setQPSubParts( pcCU->getRefQP( 0 ), 0, uiDepth ); // set QP to default QP
    }
  }
}

Void TEncCu::xCopyAMVPInfo (AMVPInfo* pSrc, AMVPInfo* pDst)
{
  pDst->iN = pSrc->iN;
  for (Int i = 0; i < pSrc->iN; i++)
  {
    pDst->m_acMvCand[i] = pSrc->m_acMvCand[i];
  }
}
Void TEncCu::xCopyYuv2Pic(TComPic* rpcPic, UInt uiCUAddr, UInt uiAbsPartIdx, UInt uiDepth, UInt uiSrcDepth, TComDataCU* pcCU, UInt uiLPelX, UInt uiTPelY )
{
  UInt uiRPelX   = uiLPelX + (g_uiMaxCUWidth>>uiDepth)  - 1;
  UInt uiBPelY   = uiTPelY + (g_uiMaxCUHeight>>uiDepth) - 1;
  TComSlice * pcSlice = pcCU->getPic()->getSlice(pcCU->getPic()->getCurrSliceIdx());
  Bool bSliceStart = pcSlice->getSliceSegmentCurStartCUAddr() > rpcPic->getPicSym()->getInverseCUOrderMap(pcCU->getAddr())*pcCU->getPic()->getNumPartInCU()+uiAbsPartIdx && 
    pcSlice->getSliceSegmentCurStartCUAddr() < rpcPic->getPicSym()->getInverseCUOrderMap(pcCU->getAddr())*pcCU->getPic()->getNumPartInCU()+uiAbsPartIdx+( pcCU->getPic()->getNumPartInCU() >> (uiDepth<<1) );
  Bool bSliceEnd   = pcSlice->getSliceSegmentCurEndCUAddr() > rpcPic->getPicSym()->getInverseCUOrderMap(pcCU->getAddr())*pcCU->getPic()->getNumPartInCU()+uiAbsPartIdx && 
    pcSlice->getSliceSegmentCurEndCUAddr() < rpcPic->getPicSym()->getInverseCUOrderMap(pcCU->getAddr())*pcCU->getPic()->getNumPartInCU()+uiAbsPartIdx+( pcCU->getPic()->getNumPartInCU() >> (uiDepth<<1) );
  if(!bSliceEnd && !bSliceStart && ( uiRPelX < pcSlice->getSPS()->getPicWidthInLumaSamples() ) && ( uiBPelY < pcSlice->getSPS()->getPicHeightInLumaSamples() ) )
  {
    UInt uiAbsPartIdxInRaster = g_auiZscanToRaster[uiAbsPartIdx];
    UInt uiSrcBlkWidth = rpcPic->getNumPartInWidth() >> (uiSrcDepth);
    UInt uiBlkWidth    = rpcPic->getNumPartInWidth() >> (uiDepth);
    UInt uiPartIdxX = ( ( uiAbsPartIdxInRaster % rpcPic->getNumPartInWidth() ) % uiSrcBlkWidth) / uiBlkWidth;
    UInt uiPartIdxY = ( ( uiAbsPartIdxInRaster / rpcPic->getNumPartInWidth() ) % uiSrcBlkWidth) / uiBlkWidth;
    UInt uiPartIdx = uiPartIdxY * ( uiSrcBlkWidth / uiBlkWidth ) + uiPartIdxX;
    m_ppcRecoYuvBest[uiSrcDepth]->copyToPicYuv( rpcPic->getPicYuvRec (), uiCUAddr, uiAbsPartIdx, uiDepth - uiSrcDepth, uiPartIdx);
  }
  else
  {
    UInt uiQNumParts = ( pcCU->getPic()->getNumPartInCU() >> (uiDepth<<1) )>>2;

    for ( UInt uiPartUnitIdx = 0; uiPartUnitIdx < 4; uiPartUnitIdx++, uiAbsPartIdx+=uiQNumParts )
    {
      UInt uiSubCULPelX   = uiLPelX + ( g_uiMaxCUWidth >>(uiDepth+1) )*( uiPartUnitIdx &  1 );
      UInt uiSubCUTPelY   = uiTPelY + ( g_uiMaxCUHeight>>(uiDepth+1) )*( uiPartUnitIdx >> 1 );

      Bool bInSlice = rpcPic->getPicSym()->getInverseCUOrderMap(pcCU->getAddr())*pcCU->getPic()->getNumPartInCU()+uiAbsPartIdx+uiQNumParts > pcSlice->getSliceSegmentCurStartCUAddr() && 
        rpcPic->getPicSym()->getInverseCUOrderMap(pcCU->getAddr())*pcCU->getPic()->getNumPartInCU()+uiAbsPartIdx < pcSlice->getSliceSegmentCurEndCUAddr();
      if(bInSlice&&( uiSubCULPelX < pcSlice->getSPS()->getPicWidthInLumaSamples() ) && ( uiSubCUTPelY < pcSlice->getSPS()->getPicHeightInLumaSamples() ) )
      {
        xCopyYuv2Pic( rpcPic, uiCUAddr, uiAbsPartIdx, uiDepth+1, uiSrcDepth, pcCU, uiSubCULPelX, uiSubCUTPelY );   // Copy Yuv data to picture Yuv
      }
    }
  }
}

Void TEncCu::xCopyYuv2Tmp( UInt uiPartUnitIdx, UInt uiNextDepth )
{
  UInt uiCurrDepth = uiNextDepth - 1;
  m_ppcRecoYuvBest[uiNextDepth]->copyToPartYuv( m_ppcRecoYuvTemp[uiCurrDepth], uiPartUnitIdx );
}

/** Function for filling the PCM buffer of a CU using its original sample array 
 * \param pcCU pointer to current CU
 * \param pcOrgYuv pointer to original sample array
 * \returns Void
 */
Void TEncCu::xFillPCMBuffer     ( TComDataCU*& pCU, TComYuv* pOrgYuv )
{

  UInt   width        = pCU->getWidth(0);
  UInt   height       = pCU->getHeight(0);

  Pel*   pSrcY = pOrgYuv->getLumaAddr(0, width); 
  Pel*   pDstY = pCU->getPCMSampleY();
  UInt   srcStride = pOrgYuv->getStride();

  for(Int y = 0; y < height; y++ )
  {
    for(Int x = 0; x < width; x++ )
    {
      pDstY[x] = pSrcY[x];
    }
    pDstY += width;
    pSrcY += srcStride;
  }

  Pel* pSrcCb       = pOrgYuv->getCbAddr();
  Pel* pSrcCr       = pOrgYuv->getCrAddr();;

  Pel* pDstCb       = pCU->getPCMSampleCb();
  Pel* pDstCr       = pCU->getPCMSampleCr();;

  UInt srcStrideC = pOrgYuv->getCStride();
  UInt heightC   = height >> 1;
  UInt widthC    = width  >> 1;

  for(Int y = 0; y < heightC; y++ )
  {
    for(Int x = 0; x < widthC; x++ )
    {
      pDstCb[x] = pSrcCb[x];
      pDstCr[x] = pSrcCr[x];
    }
    pDstCb += widthC;
    pDstCr += widthC;
    pSrcCb += srcStrideC;
    pSrcCr += srcStrideC;
  }
}

#if ADAPTIVE_QP_SELECTION
/** Collect ARL statistics from one block
  */
Int TEncCu::xTuCollectARLStats(TCoeff* rpcCoeff, Int* rpcArlCoeff, Int NumCoeffInCU, Double* cSum, UInt* numSamples )
{
  for( Int n = 0; n < NumCoeffInCU; n++ )
  {
    Int u = abs( rpcCoeff[ n ] );
    Int absc = rpcArlCoeff[ n ];

    if( u != 0 )
    {
      if( u < LEVEL_RANGE )
      {
        cSum[ u ] += ( Double )absc;
        numSamples[ u ]++;
      }
      else 
      {
        cSum[ LEVEL_RANGE ] += ( Double )absc - ( Double )( u << ARL_C_PRECISION );
        numSamples[ LEVEL_RANGE ]++;
      }
    }
  }

  return 0;
}

/** Collect ARL statistics from one LCU
 * \param pcCU
 */
Void TEncCu::xLcuCollectARLStats(TComDataCU* rpcCU )
{
  Double cSum[ LEVEL_RANGE + 1 ];     //: the sum of DCT coefficients corresponding to datatype and quantization output
  UInt numSamples[ LEVEL_RANGE + 1 ]; //: the number of coefficients corresponding to datatype and quantization output

  TCoeff* pCoeffY = rpcCU->getCoeffY();
  Int* pArlCoeffY = rpcCU->getArlCoeffY();

  UInt uiMinCUWidth = g_uiMaxCUWidth >> g_uiMaxCUDepth;
  UInt uiMinNumCoeffInCU = 1 << uiMinCUWidth;

  memset( cSum, 0, sizeof( Double )*(LEVEL_RANGE+1) );
  memset( numSamples, 0, sizeof( UInt )*(LEVEL_RANGE+1) );

  // Collect stats to cSum[][] and numSamples[][]
  for(Int i = 0; i < rpcCU->getTotalNumPart(); i ++ )
  {
    UInt uiTrIdx = rpcCU->getTransformIdx(i);

    if(rpcCU->getPredictionMode(i) == MODE_INTER)
    if( rpcCU->getCbf( i, TEXT_LUMA, uiTrIdx ) )
    {
      xTuCollectARLStats(pCoeffY, pArlCoeffY, uiMinNumCoeffInCU, cSum, numSamples);
    }//Note that only InterY is processed. QP rounding is based on InterY data only.
   
    pCoeffY  += uiMinNumCoeffInCU;
    pArlCoeffY  += uiMinNumCoeffInCU;
  }

  for(Int u=1; u<LEVEL_RANGE;u++)
  {
    m_pcTrQuant->getSliceSumC()[u] += cSum[ u ] ;
    m_pcTrQuant->getSliceNSamples()[u] += numSamples[ u ] ;
  }
  m_pcTrQuant->getSliceSumC()[LEVEL_RANGE] += cSum[ LEVEL_RANGE ] ;
  m_pcTrQuant->getSliceNSamples()[LEVEL_RANGE] += numSamples[ LEVEL_RANGE ] ;
}
#endif


//! \}
