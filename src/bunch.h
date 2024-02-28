#ifndef BUNCH_H
#define BUNCH_H

#include "wrap.h"
#include <Axis.hpp>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

class Palette{
public:
  Palette() = default;
  ~Palette() = default;
  
  Palette& SetPalette( std::vector<Style> palette ){ styles_ = std::move(palette); return *this; }
  const std::vector<Style>& GetPalette(){ return styles_; }
  template< class T >
  void PaintObjects( std::vector< Wrap<T> >& objects ) const { std::for_each( objects.begin(), objects.end(), 
    [this, i=0]( auto& obj ) mutable {
      obj.SetStyle( styles_.at(i) );
      ++i;
    } ); }
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
protected:
  std::vector<Style> styles_{};
};

template<typename T>
class Bunch {
public:
  Bunch<T>() = default;
  Bunch<T>( const Bunch<T>& other) = default;
  Bunch<T>( Bunch<T>&& other) = default;
  Bunch<T>& operator=( const Bunch<T>& other) = default;
  Bunch<T>& operator=( Bunch<T>&& ) = default;

  ~Bunch<T>() = default;

  template<typename... Args>
  Bunch& AddToBunch( std::string title, Args... args ){ 
    bunch_.emplace_back( title, args... ); 
    return *this;
  }
  Palette& GetPalette() { return palette_; }
  template<typename Func>
  Bunch& Perform( Func function ){ std::for_each( bunch_.begin(), bunch_.end(), function ); return *this; }
  std::vector<Wrap<T>>& operator*(){
    palette_.PaintObjects(bunch_);
    return bunch_; 
  }
  Wrap<T>& operator[]( size_t idx ){ return bunch_.at(idx); }
  
  TLegend* MakeLegend(std::vector<double> position = {}){
    std::vector<std::string> bunch_titles_;
    std::for_each( bunch_.begin(), bunch_.end(), [&bunch_titles_]( const Wrap<T>& obj ) mutable { 
      bunch_titles_.emplace_back( obj.GetTitle() ); 
    } );
    auto leg = palette_.MakeLegend( bunch_titles_, position );
    return leg;
  }

private:
  Palette palette_{};
  std::vector< Wrap<T> > bunch_{};
};

class DoubleDifferential {
public:
  DoubleDifferential() = default;
  template<typename ...Args>
  DoubleDifferential(Args... args) : base_correlation_{ args... }{ }
  DoubleDifferential( const DoubleDifferential& ) = default;
  DoubleDifferential& operator=( const DoubleDifferential& ) = default;
  DoubleDifferential( DoubleDifferential&& ) = default;
  DoubleDifferential& operator=( DoubleDifferential&& ) = default;
  ~DoubleDifferential() = default;
  template<typename Func>
  DoubleDifferential& Perform( const Func& function ){ 
    function( base_correlation_ ); 
    return *this; 
  }
  DoubleDifferential& SetSliceAxis(Qn::AxisD slie_axis){ slice_axis_ = std::move(slie_axis); return *this; }
  DoubleDifferential& SetProjectionAxis(Qn::AxisD projection_axis){ projection_axis_ = std::move(projection_axis); return *this; }
  DoubleDifferential& SetPalette( const std::vector<Style>& styles ){
    palette_.SetPalette(styles);
    return *this;
  }
  TLegend* MakeLegend(const std::string& slice_var, const std::string& units="", std::vector<double> position = {} ){
    auto bin_edges = slice_axis_.GetBinEdges();
    std::vector<std::string> captions{};
    for( int i = 0; i < bin_edges.size() - 1; ++i ){
      auto stream = std::ostringstream{};
      stream << std::setprecision(2) << slice_var << " " << bin_edges.at(i) << "-" << bin_edges.at(i+1) << " " << units;
      captions.emplace_back( stream.str() );
    }
    return palette_.MakeLegend( captions, position );
  }
  std::vector<Wrap<Correlation>>& operator*(){
    auto n_projections = slice_axis_.GetNBins();
    auto slice_var = slice_axis_.Name();
    auto bin_edges = slice_axis_.GetBinEdges();
    auto projection_var = projection_axis_.Name();
    base_correlation_
      .Rebin( {slice_axis_, projection_axis_} )
      .Project( {slice_var, projection_var} );
    result_correlations_.reserve(n_projections);
    for( int i=0; i<n_projections; ++i ){
      result_correlations_.emplace_back( base_correlation_ );
      
      (*result_correlations_.back()).Rebin( {{slice_var, 1, bin_edges.at(i), bin_edges.at(i+1)} } )
        .Project( {projection_var} );
      palette_.PaintObjects( result_correlations_ );
    }
    return result_correlations_;
  }
private:
  Palette palette_;
  Correlation base_correlation_;
  std::vector< Wrap<Correlation> > result_correlations_;

  Qn::AxisD slice_axis_;
  Qn::AxisD projection_axis_;
};

template<typename T>
class RatioBuilder {
public:
  template<typename... Args>
  RatioBuilder( Args... args ) : reference_(args...) {}
  ~RatioBuilder() = default;
  template<typename ...Args>
  RatioBuilder& AddToBunch( std::string title, Args... args ){ results_.emplace_back( title, args...); return *this; }
  std::vector<Wrap<T>>& GetResults(){
    palette_.PaintObjects( results_ );
    return results_;
  }
  std::vector<Wrap<T>>& GetRatios(){
    if( ratios_.empty() ){
      int i=0;
      ratios_.reserve( results_.size() );
      std::for_each( results_.begin(), results_.end(), [this, i=0](auto& obj) mutable {
        ratios_.emplace_back( obj.Divide( reference_ )  );
        ++i;
      } );
    }
    palette_.PaintObjects( ratios_ );
    return ratios_;
  }
  template<typename Func>
  RatioBuilder<T>& Perform( const Func& function ){ 
    std::for_each( results_.begin(), results_.end(), function ); 
    function( reference_ ); 
    return *this; 
  }
  RatioBuilder<T>& SetPalette( const Style& ref_style, std::vector<Style> corr_styles ){
    reference_.SetStyle(ref_style);
    palette_.SetPalette(corr_styles);
    return *this;
  }
  Wrap<T>& GetReference(){ return reference_; }
  TLegend* MakeLegend( std::vector<double> position = {} ){
    std::vector<std::string> result_titles{};
    result_titles.reserve( results_.size() );
    std::for_each( results_.begin(), results_.end(), [&result_titles] (const Wrap<T>& obj) { result_titles.push_back( obj.GetTitle() ); } );
    auto leg = palette_.MakeLegend( result_titles, position );    
    auto graph = new TGraph(0);
      graph->SetLineColor( reference_.GetStyle().color_ );
      graph->SetMarkerColor( reference_.GetStyle().color_ );
      if( reference_.GetStyle().marker_ >= 0 ){
        graph->SetMarkerStyle( reference_.GetStyle().marker_ );
        leg->AddEntry( graph, reference_.GetTitle().c_str(), "P" );
      }
      if( reference_.GetStyle().marker_ < 0 ){
        graph->SetLineStyle( reference_.GetStyle().marker_ );
        leg->AddEntry( graph, reference_.GetTitle().c_str(), "L" );
      }
      return leg;
  }
  
private:
  Palette palette_;
  Wrap<T> reference_;
  std::vector<Wrap<Correlation>> results_;
  std::vector<Wrap<Correlation>> ratios_;
};

#endif // BUNCH_H
