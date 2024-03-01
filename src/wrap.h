#ifndef WRAP_H
#define WRAP_H

#include <algorithm>
#include <cstddef>
#include <memory>
#include <queue>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <bits/fs_fwd.h>
#include <algorithm>

#include <RtypesCore.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TFile.h>

#include <DataContainer.hpp>
#include <DataContainerHelper.hpp>
#include <Axis.hpp>

struct Style{
  Style(){}
  Style& SetColor( int color ){ color_ = color; return *this; }
  Style& SetMarker( int marker ){ marker_ = marker; return *this; }
  int color_{kBlack};
  int marker_{kFullCircle};
};

template<typename T> 
class Wrap{
public:
  Wrap() = default;
  template<typename... Args>
  Wrap(std::string title, Args... args) : title_{std::move(title)}, obj_( args... ) {}
  Wrap( T obj ) : obj_( std::move(obj) ) {}
  Wrap( const Wrap& other ) : 
    title_{other.title_}, 
    points_{ nullptr },
    style_( other.style_ ),
    obj_{other.obj_}{}
  Wrap& operator=(const Wrap& other){
    title_ = other.title_; 
    points_.reset();
    style_ = other.style_;
    obj_ = other.obj_;
  }
  Wrap( Wrap&& ) = default;
  Wrap& operator=(Wrap&&) = default;
  ~Wrap() = default;
  T& operator*(){
    return obj_;
  }
  T* operator->(){
    return &obj_;
  }
  const std::string& GetTitle() const { return title_; }
  Wrap& Fit( TF1* function ){
    UpdatePoints();
    points_->Fit( function );
    auto fitted = dynamic_cast<TF1*>(points_->GetListOfFunctions()->First());
    function->SetLineColor( style_.color_ );
    fit_ = function;
    return *this;
  }
  TGraph* ReleasePoints(){ UpdatePoints(); return points_.release(); }
  TF1* GetFit() { return fit_; };
  Wrap& SetStyle( Style style ){ style_ = std::move(style); return *this; }
  void UpdatePoints(){ 
    if( !points_ )
      points_.reset( obj_.GetPoints() );
    points_->SetLineColor(style_.color_);
    points_->SetMarkerColor(style_.color_);
    if( style_.marker_ >= 0 )
      points_->SetMarkerStyle(style_.marker_);
    if( style_.marker_ < 0 )
      points_->SetLineStyle( abs(style_.marker_) );
  }
  const Style& GetStyle() const { return style_; }
  Wrap<T> Divide( const Wrap<T>& other) const {
    return Wrap<T>{ obj_.Divide( other.obj_ ) };
  }
private:
  T obj_;
  std::string title_{};
  std::unique_ptr<TGraph> points_{};
  TF1* fit_{nullptr};
  Style style_{};
};

template<class T>
class Systematics{
public:
  template<typename... Args>
  Systematics<T>( Args... args ) : object_(args...) {}
  TGraph* GetSystematics(){ if(!systamatics_) UpdatePoints(); return systamatics_.get(); }
  TGraph* ReleaseSystematics(){ if(!systamatics_) UpdatePoints(); return systamatics_.release(); }
  void UpdatePoints(){
    object_.UpdatePoints();
    systamatics_.reset( object_->GetSystematics() );
    systamatics_->SetLineColor( object_.GetStyle().color_ );
    if( object_.GetStyle().marker_ >= 0 ){
      systamatics_->SetMarkerStyle( object_.GetStyle().marker_ );
      systamatics_->SetMarkerColor( object_.GetStyle().color_ );
    } else {
      systamatics_->SetLineStyle( abs(object_.GetStyle().marker_) );
    }
    systamatics_->SetFillColorAlpha(object_.GetStyle().color_, 0.1);
  }
  Wrap<T>* operator->(){
    return &object_;
  }

private:
  Wrap<T> object_;
  std::unique_ptr<TGraph> systamatics_;
};

class Correlation{
public:
  Correlation() = default;
  Correlation( 
    std::string str_file_name, 
    std::vector<std::string> vec_objects,
    std::vector<double> vec_weihts = {}
  ){
    auto file = std::make_unique<TFile>( str_file_name.c_str(), "READ" );
    Qn::DataContainerStatCalculate* ptr{nullptr};
    std::queue<Qn::DataContainerStatCalculate> correlations;
    int i=0;
    for( auto name : vec_objects ){
      auto weight = 1.;
      if( !vec_weihts.empty() )
        weight = vec_weihts.at(i);
      file->GetObject( name.c_str(), ptr );
      if( !ptr ){
        Qn::DataContainerStatCalculate* ptr_stat_collect{nullptr};
        file->GetObject( name.c_str(), ptr_stat_collect );
        if( !ptr_stat_collect )
            throw std::runtime_error( std::string("No object in a file").append( " " ).append( name ) );
        correlations.push( Qn::DataContainerStatCalculate(*ptr_stat_collect)*weight );
      }
      correlations.push( Qn::DataContainerStatCalculate(*ptr)*weight );
      ++i;
    }
    average_ = correlations.front();
    averaging_objects_.emplace_back( average_ );
    correlations.pop();
    auto* list_merge = new TList;
    while( !correlations.empty() ) {
      auto* to_merge = new Qn::DataContainerStatCalculate( correlations.front() );
      list_merge->Add(to_merge);
      averaging_objects_.emplace_back( correlations.front() );
      correlations.pop();
    }
    average_.Merge(list_merge);
    
  }
  Correlation( Qn::DataContainerStatCalculate corr ) : average_{ std::move(corr) } {}
  Correlation(const Correlation& ) = default;
  Correlation& operator=(const Correlation& ) = default;
  Correlation(Correlation&& ) = default;
  Correlation& operator=(Correlation&& ) = default;
  ~Correlation() = default;
  TGraphErrors* GetPoints() {
    average_.SetErrors(Qn::Stat::ErrorType::BOOTSTRAP);
    auto graph = Qn::ToTGraph( average_ );
    return graph;
  }
  TGraphErrors* GetSystematics(){
    auto graph = Qn::ToTGraph( average_ );
    std::vector<Qn::DataContainerStatCalculate> variations;
    for( const auto& cont : averaging_objects_ ){
      variations.emplace_back( average_ - cont );
    }
    for( size_t i=0; i<graph->GetN(); ++i ){
      auto x_hi = average_.GetAxes().front().GetUpperBinEdge(i);
      auto x_lo = average_.GetAxes().front().GetLowerBinEdge(i);
      auto x_err = fabs( x_hi - x_lo ) / 4;
      std::vector<double> vec_variations;
      auto y_err = fabs(variations.front().At(i).Mean());
      for( const auto& cont : variations ){
        auto err = fabs( cont.At(i).Mean() );
        if( y_err < err )
          y_err=err;
      }
      graph->SetPointError( i, x_err, y_err );
    }
    return graph;
  }
  Correlation& Scale(double scale){
    average_ = average_ * scale;
    for( auto& corr : averaging_objects_ ){
      corr = corr*scale;
    }
    return *this;
  }
  Correlation& Rebin( const std::vector<Qn::AxisD>& axes ){
    for( const auto& axis : axes ){
      average_ = average_.Rebin(axis);
      for( auto& corr : averaging_objects_ ){
        corr = corr.Rebin( axis );
      }
    }
    return *this;
  }
  Correlation& Project( const std::vector<std::string>& axes ){
    average_ = average_.Projection( axes );
    for( auto& corr : averaging_objects_ ){
        corr = corr.Projection( axes );
    }
    return *this;
  }
  Correlation Divide( const Correlation& other ) const {
    auto result_correlation = average_ / other.average_;
    return Correlation(result_correlation);
  }
private:
  Qn::DataContainerStatCalculate  average_{};
  std::vector<Qn::DataContainerStatCalculate>  averaging_objects_{};
};

class Histogram{
public:
  Histogram( 
    const std::string& str_file_name, 
    const std::vector<std::string>& vec_objects )
  {
    auto file = std::make_unique<TFile>( str_file_name.c_str(), "READ" );
    std::queue<TH1*> histograms;
    TH1* ptr;
    for( const auto& name : vec_objects ){
      file->GetObject( name.c_str(), ptr );
      if( !ptr )
        throw std::runtime_error( std::string("No object in a file").append( " " ).append( name ) );
      histograms.push( ptr );
    }
    auto new_name = std::string(histograms.front()->GetName()) + "_copy";
    histogram_ = std::unique_ptr<TH1>( dynamic_cast<TH1*>(histograms.front()->Clone( new_name.c_str() )) );
    histograms.pop();
    auto* list_merge = new TList;
    while( !histograms.empty() ) {
      auto* to_merge = histograms.front();
      list_merge->Add(to_merge);
      histograms.pop();
    }
    histogram_->Merge( list_merge );
  }
  Histogram(const Histogram& other) : 
    histogram_{ dynamic_cast<TH1*>(other.histogram_->Clone( std::data( std::string{ histogram_->GetName() }.append("_copy") ) )) }{};
  Histogram& operator=(const Histogram& other){
    histogram_.reset( dynamic_cast<TH1*>( other.histogram_->Clone( std::data( std::string{ histogram_->GetName() }.append("_copy") ) )) );
    return *this;
  };
  Histogram(Histogram&& ) = default;
  Histogram& operator=(Histogram&& ) = default;
  Histogram( TH1* hist ) : histogram_( std::unique_ptr<TH1>( hist ) ) {}
  TGraphErrors* GetPoints() const {
    auto n_bins = histogram_->GetNbinsX();
    auto x_axis = std::vector<double>{};
    auto y_axis = std::vector<double>{};
    auto y_error = std::vector<double>{};
    for( int i=0; i<n_bins; ++i ){
      x_axis.push_back( histogram_->GetBinCenter(i+1) );
      y_axis.push_back( histogram_->GetBinContent(i+1) );
      y_error.push_back( histogram_->GetBinError(i+1) );
    }
    auto graph = new TGraphErrors(n_bins, x_axis.data(), y_axis.data(), nullptr, y_axis.data() );
    return graph;
  }
  Histogram& Scale(double scale){
    histogram_->Scale(scale);
    return *this;
  }
  Histogram& Rebin( size_t n_group ){
    histogram_->RebinX(n_group);
    return *this;
  }
  Histogram Divide( const Histogram& other ) const {
    auto result_histogram = dynamic_cast<TH1*>(histogram_->Clone(std::data(std::string("divide").append( histogram_->GetName() )) ));
    result_histogram->Divide( other.histogram_.get() );
    return Histogram(result_histogram);
  }
private:
  std::unique_ptr<TH1> histogram_{};
};

class Graph{
  Graph( std::string str_file_name, std::string str_obj ){
    auto file = std::make_unique<TFile>( str_file_name.c_str(), "READ" );
    TGraph* ptr{nullptr};
    file->GetObject(str_obj.c_str(), ptr);
    graph_ = std::unique_ptr<TGraphErrors>( dynamic_cast<TGraphErrors*>(ptr->Clone()) );
  }
  Graph( const Histogram&  histogram ) : graph_{ std::unique_ptr<TGraphErrors>( dynamic_cast<TGraphErrors*>( histogram.GetPoints() ) )}{ }
  Graph( Correlation&  correlation ) : graph_ {std::unique_ptr<TGraphErrors>( dynamic_cast<TGraphErrors*>( correlation.GetPoints() ) )}{ }
  Graph& Scale(double scale){
    for( size_t i=0; i<graph_->GetN(); ++i ){
      auto x = graph_->GetPointX(i);
      auto y = graph_->GetPointY(i); 
      auto y_err = graph_->GetErrorY(i); 
      graph_->SetPoint( i, x,  y*scale );
      graph_->SetPointError( i, 0,  y_err*scale );
    }
    return *this;
  }
  Graph& ScaleXaxis(double scale){
    for( size_t i=0; i<graph_->GetN(); ++i ){
      auto x = graph_->GetPointX(i);
      auto y = graph_->GetPointY(i); 
      auto y_err = graph_->GetErrorY(i); 
      graph_->SetPoint( i, x*scale, y );
      graph_->SetPointError( i, 0, y_err );
    }
    return *this;
  }
  Graph& SetXAxis( std::vector<double> x_axis ){
    int i=0;
    for( auto x : x_axis ){
      graph_->SetPointX( i, x);
      ++i;
    }
    return *this;
  }
private:
  std::unique_ptr<TGraphErrors> graph_;
};

#endif // WRAP_H