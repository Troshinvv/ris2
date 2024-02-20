#ifndef BUNCH_H
#define BUNCH_H

#include "wrap.h"
#include <Axis.hpp>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <algorithm>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

template<typename T, typename Base>
class Palette{
public:
  Palette() = default;
  ~Palette() = default;
  template<typename Func>
  Base& Perform( Func function ){ std::for_each( bunch_.begin(), bunch_.end(), function ); return static_cast<Base&>(*this); }
  Base& SetStyles( std::vector<Style> styles ){ styles_ = std::move( styles ); return static_cast<Base&>(*this); }
  TLegend* MakeLegend( const std::vector<std::string>& captions, std::vector<double> position = {} ){
    int i=0;
    TLegend* leg{};
    if( position.size() == 4 )
      leg = new TLegend{ position.at(0), position.at(1), position.at(2), position.at(3) };
    else{
      leg = new TLegend{};
    }
    for( auto& titile : captions ){
      auto graph = new TGraph(0);
      graph->SetLineColor( styles_.at(i).color_ );
      graph->SetMarkerColor( styles_.at(i).color_ );
      if( styles_.at(i).marker_ >= 0 ){
        graph->SetMarkerStyle( styles_.at(i).marker_ );
        leg->AddEntry( graph, titile.c_str(), "P" );
      }
      if( styles_.at(i).marker_ < 0 ){
        graph->SetLineStyle( styles_.at(i).marker_ );
        leg->AddEntry( graph, titile.c_str(), "L" );
      }
      ++i;
    }
    return leg;
  }
  std::vector<Wrap<T>>& operator*(){
    std::for_each( bunch_.begin(), bunch_.end(), [this, idx=0]( Wrap<T>& obj ) mutable { obj.SetStyle( styles_.at(idx) ); idx++; });
    return bunch_; 
  }
  T* operator[]( size_t idx ){ return bunch_.at(idx); }
protected:
  std::vector< Wrap<T> > bunch_{};
  std::vector<Style> styles_{};
};

template<typename T>
class Bunch : public Palette<T, Bunch<T>>{
public:
  Bunch() = default;
  ~Bunch() = default;
  template<typename... Args>
  Bunch& AddToBunch( std::string title, Args... args ){ 
    Palette<T, Bunch<T> >::bunch_.emplace_back( title, args... ); return *this; 
    return *this;
  }
};

class DoubleDifferential : public Palette<Correlation, DoubleDifferential> {
public:
  DoubleDifferential() = default;
  template<typename ...Args>
  DoubleDifferential(Args... args) : base_correlation_{ args... }{ }
  DoubleDifferential( const DoubleDifferential& ) = default;
  DoubleDifferential& operator=( const DoubleDifferential& ) = default;
  DoubleDifferential( DoubleDifferential&& ) = default;
  DoubleDifferential& operator=( DoubleDifferential&& ) = default;
  ~DoubleDifferential() = default;
  DoubleDifferential& Rebin( const std::vector<Qn::AxisD>& axes ){ base_correlation_.Rebin(axes); return *this; }
  DoubleDifferential& Project( const std::vector<std::string>& axes ){ base_correlation_.Project(axes); return *this; }
  DoubleDifferential& SetSliceAxis(Qn::AxisD slie_axis){ slice_axis_ = std::move(slie_axis); return *this; }
  DoubleDifferential& SetProjectionAxis(Qn::AxisD projection_axis){ projection_axis_ = std::move(projection_axis); return *this; }
  TLegend* MakeLegend(const std::string& slice_var, const std::string& units="", std::vector<double> position = {} ){
    auto bin_edges = slice_axis_.GetBinEdges();
    std::vector<std::string> captions{};
    captions.reserve( bin_edges.size() - 1 );
    for( int i=0; i<bin_edges.size()-1; ++i ){
      auto stream = std::ostringstream{};
      stream << std::setprecision(2) << bin_edges.at(i) << "-" << bin_edges.at(i+1) << " " << units;
      captions.emplace_back( stream.str() );
    }
    return Palette<Correlation, DoubleDifferential>::MakeLegend(captions, position);
  }
  std::vector<Wrap<Correlation>>& operator*(){
    auto n_projections = slice_axis_.GetNBins();
    auto slice_var = slice_axis_.Name();
    auto bin_edges = slice_axis_.GetBinEdges();
    auto projection_var = projection_axis_.Name();
    base_correlation_
      .Rebin( {slice_axis_, projection_axis_} )
      .Project( {slice_var, projection_var} );
    bunch_.reserve(n_projections);
    for( int i=0; i<n_projections-1; ++i ){
      bunch_.push_back( base_correlation_ );
      bunch_.back()
        .Rebin({{slice_var, 1, bin_edges.at(i), bin_edges.at(i+1)} } )
        .Project( {projection_var} );
    }
    return Palette<Correlation, DoubleDifferential>::operator*();
  }
private:
  Correlation base_correlation_;
  Qn::AxisD slice_axis_;
  Qn::AxisD projection_axis_;
};

template<typename T>
class RatioBuilder : public Palette<T, RatioBuilder<T>>{
public:
  template<typename... Args>
  RatioBuilder( Args... args ) : reference_(args...) {}
  ~RatioBuilder() = default;
  template<typename ...Args>
  RatioBuilder& AddToBunch( std::string title, Args... args ){ 
    Palette<T, RatioBuilder<T> >::bunch_.emplace_back( title, args... ); return *this; 
    return *this;
  }
  std::vector<Wrap<T>>& GetResults(){
    return Palette<T, RatioBuilder<T>>::operator*();
  }
  std::vector<Wrap<T>>& GetRatios(){
    if( ratios_.empty() ){
      int i=0;
      ratios_.reserve( Palette<T, RatioBuilder<T>>::bunch_.size() );
      std::for_each( Palette<T, RatioBuilder<T>>::bunch_.begin(), Palette<T, RatioBuilder<T>>::bunch_.end(), [this, i=0](auto& obj) mutable {
        ratios_.emplace_back( obj.Divide(reference_) );
        ratios_.back().SetStyle(Palette<T, RatioBuilder<T>>::styles_.at(i));
        ++i;
      } );
    }
    return ratios_;
  }
  template<typename Func>
  RatioBuilder<T>& Perform( Func function ){ Palette<T, RatioBuilder<T>>::Perform(function); function( reference_ ); return *this; }
  Wrap<T>& GetReference(){ return reference_; }
  TLegend* MakeLegend( const std::string ref_name, const std::vector<std::string>& captions, std::vector<double> position = {} ){
    auto leg = Palette<T, RatioBuilder<T>>::MakeLegend( captions, position );
    auto graph = new TGraph(0);
    graph->SetLineColor( reference_.GetStyle().color_ );
    graph->SetMarkerColor( reference_.GetStyle().color_ );
    if( reference_.GetStyle().marker_ >= 0 ){
      graph->SetMarkerStyle( reference_.GetStyle().marker_ );
      leg->AddEntry( graph, ref_name.c_str(), "P" );
    }
    if( reference_.GetStyle().marker_ < 0 ){
      graph->SetLineStyle( reference_.GetStyle().marker_ );
      leg->AddEntry( graph, ref_name.c_str(), "L" );
    };
    return leg;
  }
private:
  Wrap<T> reference_;
  std::vector< Wrap<T> > ratios_;
};

#endif // BUNCH_H
